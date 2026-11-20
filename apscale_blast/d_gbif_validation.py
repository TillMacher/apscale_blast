from pygbif import occurrences as occ
from pygbif import species as pygbif_species
import math
import operator
import re
import time
from pathlib import Path
import pandas as pd
from tqdm import tqdm
import pycountry
import math
import plotly.graph_objects as go
import logging

logging.getLogger("kaleido").setLevel(logging.WARNING)
logging.getLogger("choreographer").setLevel(logging.WARNING)

## Flags added by the GBIF validation
FLAG_GBIF_MOST_ABUNDANT = 'F5 (GBIF: most abundant of ambiguous taxa)'
FLAG_GBIF_NO_OCCURRENCE = 'F6 (GBIF: no occurrence in region)'
FLAG_GBIF_GENUS_REPLACEMENT = 'F7 (GBIF: replaced by most abundant species of genus)'
FLAG_GBIF_SYNONYM = 'F8 (GBIF: synonym replaced by accepted name)'

def country_validation(countries):
    for code in countries:
        res = pycountry.countries.get(alpha_2=code.strip().upper())
        if not res:
            return code
    return True

def continent_validation(continent):
    allowed_countries = ['africa', 'antarctica', 'asia', 'europe', 'north_america', 'oceania', 'south_america']
    if continent not in allowed_countries:
        return continent
    else:
        return True

def lat_lon_radius_validation(lat_lon_radius):
    try:
        lat = float(lat_lon_radius.split(',')[0])
        lon = float(lat_lon_radius.split(',')[1])
        radius = int(lat_lon_radius.split(',')[2])
        return [lat, lon, radius]
    except ValueError:
        return False

def circle_wkt(lat, lon, radius_km, n_points=8):
    """WKT polygon approximating a circle of radius_km around (lat, lon), counter-clockwise for GBIF."""
    R = 6371.0088  # mean Earth radius (km)
    d = radius_km / R
    lat1, lon1 = math.radians(lat), math.radians(lon)
    coords = []
    for i in range(n_points):
        brg = -2 * math.pi * i / n_points  # decreasing bearing = counter-clockwise
        lat2 = math.asin(math.sin(lat1) * math.cos(d) + math.cos(lat1) * math.sin(d) * math.cos(brg))
        lon2 = lon1 + math.atan2(math.sin(brg) * math.sin(d) * math.cos(lat1),
                                 math.cos(d) - math.sin(lat1) * math.sin(lat2))
        coords.append(f'{math.degrees(lon2):.5f} {math.degrees(lat2):.5f}')
    coords.append(coords[0])  # close the ring
    return f'POLYGON(({", ".join(coords)}))'

def taxa_to_keys(all_species):
    # Exact match against the GBIF Backbone (same taxonomy as the occurrence index)
    # Synonyms resolve to their accepted species: all names of one species share one key
    # Names with special characters (e.g. 'Cottus asper/rhenanus', 'Salmo sp.') are skipped
    COL_XR = '7ddf754f-d193-4cc9-b351-99906754a03b'
    gbif_key_dict = {}       # species or genus name -> GBIF key
    accepted_name_dict = {}  # species name -> accepted species name
    for species in tqdm(all_species, desc='Fetching taxa keys'):
        species_name = species
        genus_name = species.split(' ')[0]

        # Only plain names (letters and spaces)
        if not re.fullmatch(r'[A-Za-z ]+', species_name):
            continue

        try:
            res = pygbif_species.name_backbone(scientificName=species_name, strict=True, checklistKey=COL_XR)
        except Exception:
            continue
        time.sleep(0.2)

        # Only accept exact matches (no fuzzy or higher-rank matches)
        match_type = res.get('diagnostics', {}).get('matchType', res.get('matchType', ''))
        if match_type != 'EXACT':
            continue

        # Accepted species (synonyms resolve to their accepted species) and its genus
        classification = {c.get('rank'): c for c in res.get('classification', [])}
        accepted = res.get('acceptedUsage') or res.get('usage') or {}
        if 'SPECIES' in classification:
            species_key = classification['SPECIES'].get('key')
            accepted_name = classification['SPECIES'].get('name', '')
        elif accepted.get('rank') == 'SPECIES':
            species_key = accepted.get('key')
            accepted_name = accepted.get('canonicalName', '')
        else:
            species_key = res.get('speciesKey')
            accepted_name = res.get('species', '')
        genus_key = classification.get('GENUS', {}).get('key') or res.get('genusKey')

        if species_key:
            accepted_name = accepted_name if accepted_name != '' else species_name
            accepted_name_dict[species_name] = accepted_name
            gbif_key_dict[species_name] = str(species_key)
            gbif_key_dict[accepted_name] = str(species_key)
            if accepted_name != species_name:
                tqdm.write(f'Synonym: {species_name} -> {accepted_name}')
        if genus_key:
            # Store under the original and the accepted genus name (synonyms can move genera)
            gbif_key_dict[genus_name] = str(genus_key)
            gbif_key_dict[accepted_name_dict.get(species_name, species_name).split(' ')[0]] = str(genus_key)
    return gbif_key_dict, accepted_name_dict

def taxon_counts(keys, batch_size=50, **region):
    COL_XR = '7ddf754f-d193-4cc9-b351-99906754a03b'
    keys = [str(k) for k in keys]
    counts = {k: 0 for k in keys}
    for i in range(0, len(keys), batch_size):
        batch = keys[i:i + batch_size]
        res = occ.search(
            taxonKey=batch,
            checklistKey=COL_XR,
            hasGeospatialIssue=False,
            facet='speciesKey',
            facetLimit=len(batch) * 5,
            limit=0,
            **region,
        )
        for c in res['facets'][0]['counts'] if res['facets'] else []:
            if c['name'] in counts:
                counts[c['name']] = c['count']
        time.sleep(1)
    return counts

def top_species_in_genus(genus_key, **region):
    COL_XR = '7ddf754f-d193-4cc9-b351-99906754a03b'
    result = {}
    res = occ.search(
        taxonKey=str(genus_key),
        checklistKey=COL_XR,
        hasGeospatialIssue=False,
        facet='speciesKey',
        facetLimit=2,          # facets are sorted by count, so these are the top two
        limit=0,
        **region,
    )
    top = res['facets'][0]['counts'] if res['facets'] else []
    if not top or (len(top) == 2 and top[0]['count'] == top[1]['count']):
        result[genus_key] = False
    else:
        match = pygbif_species.name_backbone(scientificName='', usageKey=str(top[0]['name']), checklistKey=COL_XR)
        result[genus_key] = match.get('usage', {}).get('canonicalName', False)
    time.sleep(1)
    return result

def gbif_habitat_by_key(species_key):
    HABITATS = ('freshwater', 'marine', 'terrestrial')
    res = pygbif_species.name_usage(key=species_key, data='speciesProfiles', limit=100)
    found = set()
    for profile in res.get('results', []):
        text = (profile.get('habitat') or '').lower()
        for h in HABITATS:
            if h in text or profile.get(h):
                found.add(h)
    return ', '.join(h for h in HABITATS if h in found)

def split_ambiguous_taxa(ambiguous_taxa):
    # 'Species a (99.1), Species b (98.7)' -> ['Species a', 'Species b']
    return [x.strip() for x in re.sub(r'\s*\([^)]*\)', '', ambiguous_taxa).split(',') if x.strip() != '']

def plot_map_latlonrad(lat, lon, radius_km, output_file):
    # Calculate points around the circle (1 degree latitude ≈ 111 km)
    circle_lats = []
    circle_lons = []
    for angle in range(0, 361, 5):
        a = math.radians(angle)
        circle_lats.append(lat + (radius_km / 111) * math.sin(a))
        circle_lons.append(lon + (radius_km / (111 * math.cos(math.radians(lat)))) * math.cos(a))

    # Circle
    fig = go.Figure()
    fig.add_trace(go.Scattermap(
        fill="toself",
        lon=circle_lons, lat=circle_lats,
        mode="lines",
        line={'color': "red"},
        fillcolor="rgba(255, 0, 0, 0.2)"))

    # Centre point
    fig.add_trace(go.Scattermap(
        lon=[lon], lat=[lat],
        mode="markers",
        marker={'size': 10, 'color': "red"}))

    # Rough zoom so the circle fits (bigger radius = smaller zoom)
    zoom = 12.5 - math.log2(radius_km)

    fig.update_layout(
        map={
            'style': "open-street-map",
            'center': {'lon': lon, 'lat': lat},
            'zoom': zoom},
        showlegend=False,
        margin=dict(l=0, r=0, t=0, b=0))

    fig.write_image(output_file)

def plot_map_countries(country_codes, output_file):
    iso3_codes = [pycountry.countries.get(alpha_2=code).alpha_3 for code in country_codes]

    # Fill the selected countries in red
    fig = go.Figure(go.Choropleth(
        locations=iso3_codes,
        locationmode="ISO-3",
        z=[1] * len(iso3_codes),
        colorscale=[[0, "red"], [1, "red"]],
        marker_opacity=0.6,
        marker_line_color="white",
        showscale=False))

    fig.update_layout(
        geo={
            'fitbounds': "locations",   # zoom automatically to the selected countries
            'resolution': 110,
            'showcountries': True,
            'countrycolor': "grey",
            'showland': True,
            'landcolor': "lightgrey"},
        margin=dict(l=0, r=0, t=0, b=0))

    fig.write_image(output_file)

def main(taxonomy_table_xlsx, continent, countries, lat_lon_radius, similarity_threshold, habitat_lookup, flag_7, synonyms_only=False):
    # Prepare variables
    taxonomy_table_xlsx = Path(taxonomy_table_xlsx)
    species_row = []
    flag_row = []
    occurrence_row = []
    genus_top_species_dict = {}

    # Validate user input and define the search region (not needed when only checking synonyms)
    if synonyms_only == True:
        print(' Note: Only checking for synonyms (no occurrence search).')

    elif continent != False:
        continent = continent.lower()
        res = continent_validation(continent)
        if res != True:
            print(f'Continent is not valid: {res}')
            return
        print(f' Note: Searching on continent: {continent}')
        region = {'continent': continent.upper()}

    elif countries != False:
        countries = [i.strip().upper() for i in countries.split(',')]
        res = country_validation(countries)
        if res != True:
            print(f'Country is not valid: {res}')
            return
        print(f' Note: Searching for countries: {countries}')
        region = {'country': countries}

    elif lat_lon_radius != False:
        llr = lat_lon_radius_validation(lat_lon_radius)
        if llr == False:
            print(f'Lat lon radius is not valid: {lat_lon_radius}')
            return
        print(f' Note: Searching in polygon: {llr[0]}, {llr[1]}, {llr[2]}km')
        lat, lon, radius_km = llr
        region = {'geometry': circle_wkt(lat, lon, radius_km, n_points=8)}

    else:
        print('At least one of "continent", "countries", or "latlonrad" is required')
        return

    # Read taxonomy table
    df = pd.read_excel(taxonomy_table_xlsx).fillna('')

    # Collect all species (assigned + ambiguous taxa)
    all_species = [i for i in df['Species'].unique().tolist() if i != '']
    for ambiguous_taxa in df['Ambiguous taxa'].values.tolist():
        if ambiguous_taxa != '':
            all_species.extend(split_ambiguous_taxa(ambiguous_taxa))
    all_species = sorted(set(all_species))

    # Fetch GBIF keys (species and their genera) and accepted names
    gbif_key_dict, accepted_name_dict = taxa_to_keys(all_species)

    # Replace synonyms by accepted names (keep the original assignment)
    df['Species (original)'] = df['Species']
    df['Species'] = [accepted_name_dict.get(i, i) for i in df['Species']]
    synonym_mask = df['Species'] != df['Species (original)']
    df.loc[synonym_mask, 'Genus'] = df.loc[synonym_mask, 'Species'].str.split().str[0]
    all_species = sorted({accepted_name_dict.get(i, i) for i in all_species})

    # Fetch habitats (species only)
    habitat_dict = {}
    if habitat_lookup == True:
        for species in tqdm(all_species, desc='Fetching habitats'):
            species_key = gbif_key_dict.get(species, '')
            if species_key != '':
                habitat_dict[species] = gbif_habitat_by_key(species_key)
            else:
                habitat_dict[species] = ''

    # Synonyms only: replace synonyms by accepted names and write the table (no occurrence search)
    if synonyms_only == True:
        species_row = []
        flag_row = []
        for i, row in df.iterrows():
            species = row['Species']
            flag = row['Flag']

            # Ambiguous taxa that are all synonyms of one species: assign this species (F8)
            if flag != '' and row['Similarity'] >= similarity_threshold:
                ambiguous_species = split_ambiguous_taxa(row['Ambiguous taxa'])
                ambiguous_species = list(dict.fromkeys([accepted_name_dict.get(name, name) for name in ambiguous_species]))
                if len(ambiguous_species) == 1:
                    species = ambiguous_species[0]
                    flag = f'{flag}, {FLAG_GBIF_SYNONYM}'

            # Species replaced by its accepted name (F8)
            if synonym_mask[i] == True and FLAG_GBIF_SYNONYM not in flag:
                flag = f'{flag}, {FLAG_GBIF_SYNONYM}' if flag != '' else FLAG_GBIF_SYNONYM

            species_row.append(species)
            flag_row.append(flag)

        # Write the table
        df2 = df.copy()
        df2['Species'] = species_row
        df2['Flag'] = flag_row
        mask = (df2['Genus'] == '') & (df2['Species'] != '')
        df2.loc[mask, 'Genus'] = df2.loc[mask, 'Species'].str.split().str[0]
        if habitat_lookup == True:
            df2['Habitat'] = df2['Species'].map(habitat_dict).fillna('')
        n_synonyms = sum([FLAG_GBIF_SYNONYM in i for i in flag_row])
        print(f' Note: {n_synonyms} assignments were replaced by their accepted name.')
        try:
            taxonomy_table_synonyms_xlsx = taxonomy_table_xlsx.parent / taxonomy_table_xlsx.name.replace('.xlsx', '_gbif_synonyms.xlsx')
            df2.to_excel(taxonomy_table_synonyms_xlsx, sheet_name='Taxonomy table', index=False)
        except Exception:
            taxonomy_table_synonyms_csv = taxonomy_table_xlsx.parent / taxonomy_table_xlsx.name.replace('.xlsx', '_gbif_synonyms.csv')
            df2.to_csv(taxonomy_table_synonyms_csv, index=False)
        return

    # Fetch occurrences for all species keys (batched)
    all_species_keys = sorted({gbif_key_dict[i] for i in all_species if gbif_key_dict.get(i, '') != ''})
    occurrence_dict = taxon_counts(all_species_keys, **region)

    def species_occurrences(species):
        # Number of GBIF occurrences in the region (0 if the species has no GBIF key)
        species_key = gbif_key_dict.get(species, '')
        if species_key == '':
            return 0
        return occurrence_dict.get(str(species_key), 0)

    def genus_top_species(genus):
        # Most abundant species of the genus in the region (False if none or tie); one request per genus
        if genus not in genus_top_species_dict:
            genus_key = gbif_key_dict.get(genus, '')
            if genus == '' or genus_key == '':
                genus_top_species_dict[genus] = False
            else:
                genus_top_species_dict[genus] = top_species_in_genus(genus_key, **region)[genus_key]
        return genus_top_species_dict[genus]

    # Process taxa
    for idx, row in tqdm(df.iterrows(), total=len(df), desc='Processing taxonomy table'):
        similarity = row['Similarity']
        species = row['Species']
        flag = row['Flag']

        # Only validate sequences with high similarity
        if similarity < similarity_threshold:
            species_row.append(species)
            flag_row.append(flag)
            occurrence_row.append('')
            continue

        # A) Flagged assignments: choose between the ambiguous taxa
        if flag != '':
            ambiguous_species = split_ambiguous_taxa(row['Ambiguous taxa'])
            if ambiguous_species == []:
                species_row.append(species)
                flag_row.append(flag)
                occurrence_row.append('')
                continue

            # Merge synonyms: replace by accepted names and drop duplicates (order is kept)
            ambiguous_species = list(dict.fromkeys([accepted_name_dict.get(i, i) for i in ambiguous_species]))

            # All ambiguous taxa are synonyms of one species: F8
            if len(ambiguous_species) == 1:
                species_row.append(ambiguous_species[0])
                flag_row.append(f'{flag}, {FLAG_GBIF_SYNONYM}')
                occurrence_row.append(f'{ambiguous_species[0]} (o={species_occurrences(ambiguous_species[0])})')
                continue

            # Check how often each ambiguous taxon occurs
            ambiguous_occurrence_dict = {i: species_occurrences(i) for i in ambiguous_species}

            # Rank by occurrences (second place is -1 if there is only one taxon)
            ambiguous_occurrence_sorted = sorted(ambiguous_occurrence_dict.items(), key=operator.itemgetter(1), reverse=True)
            most_occurrence = ambiguous_occurrence_sorted[0]
            second_occurrence = ambiguous_occurrence_sorted[1] if len(ambiguous_occurrence_sorted) > 1 else ('', -1)

            # No occurrences at all: F6, or F7 if enabled and the genus has one most abundant species
            if most_occurrence[1] == 0:
                top_species = genus_top_species(row['Genus']) if flag_7 == True else False
                if top_species != False:
                    species_row.append(top_species)
                    flag_row.append(f'{flag}, {FLAG_GBIF_GENUS_REPLACEMENT}')
                else:
                    species_row.append('')
                    flag_row.append(f'{flag}, {FLAG_GBIF_NO_OCCURRENCE}')

            # One clearly most abundant taxon: F5
            elif most_occurrence[1] > second_occurrence[1]:
                species_row.append(most_occurrence[0])
                flag_row.append(f'{flag}, {FLAG_GBIF_MOST_ABUNDANT}')

            # Tie: keep the original assignment
            else:
                species_row.append(species)
                flag_row.append(flag)

            occurrence_row.append(', '.join([f'{i} (o={ambiguous_occurrence_dict[i]})' for i in ambiguous_species]))
            continue

        # B) Unflagged species assignment
        if species != '':
            n_occurrences = species_occurrences(species)

            # Species occurs in the region: keep
            if n_occurrences != 0:
                species_row.append(species)
                flag_row.append(flag)

            # No occurrences: F6, or F7 if enabled and the genus has one most abundant species
            else:
                top_species = genus_top_species(row['Genus']) if flag_7 == True else False
                if top_species != False:
                    species_row.append(top_species)
                    flag_row.append(FLAG_GBIF_GENUS_REPLACEMENT)
                else:
                    species_row.append('')
                    flag_row.append(FLAG_GBIF_NO_OCCURRENCE)

            occurrence_row.append(f'{species} (o={n_occurrences})')

        # C) No species assignment: nothing to validate
        else:
            species_row.append(species)
            flag_row.append(flag)
            occurrence_row.append('')

    # Synonym flag for all rows whose species was replaced by its accepted name
    for i, replaced in enumerate(synonym_mask.tolist()):
        if replaced == True and FLAG_GBIF_SYNONYM not in flag_row[i]:
            flag_row[i] = f'{flag_row[i]}, {FLAG_GBIF_SYNONYM}' if flag_row[i] != '' else FLAG_GBIF_SYNONYM

    # Write the validated table
    df2 = df.copy()
    df2['Species'] = species_row
    df2['GBIF occ.'] = occurrence_row
    df2['Flag'] = flag_row
    mask = (df2['Genus'] == '') & (df2['Species'] != '')
    df2.loc[mask, 'Genus'] = df2.loc[mask, 'Species'].str.split().str[0]
    if habitat_lookup == True:
        df2['Habitat'] = df2['Species'].map(habitat_dict).fillna('')
    try:
        taxonomy_table_gbif_xlsx = taxonomy_table_xlsx.parent / taxonomy_table_xlsx.name.replace('.xlsx', '_gbif.xlsx')
        df2.to_excel(taxonomy_table_gbif_xlsx, sheet_name='Taxonomy table', index=False)
    except Exception:
        taxonomy_table_gbif_csv = taxonomy_table_xlsx.parent / taxonomy_table_xlsx.name.replace('.xlsx', '_gbif.csv')
        df2.to_csv(taxonomy_table_gbif_csv, index=False)

    # Plot maps
    map_file = taxonomy_table_xlsx.parent / taxonomy_table_xlsx.name.replace('.xlsx', '_map.png')
    if countries != False:
        plot_map_countries(countries, map_file)
    if lat_lon_radius != False:
        plot_map_latlonrad(lat, lon, radius_km, map_file)

if __name__ == '__main__':
    main()