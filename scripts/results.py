
# phageIDs
cities_df = pd.read_csv('cities.csv', header=None)
prophages_df = pd.read_csv('prophages.txt', sep='\t')

cities = list(map(str.upper, list(cities_df.iloc[:, 0])))
used = []

def get_city(cities_list, used):
    city = random.choice(cities)
    if city not in used:
        return city
    else:
        return get_city(cities_list, used)

phageIDs = []
for i in range(len(prophages_df)):
    city = get_city(cities, used)
    one, two, three, four = str(random.randint(0,9)), str(random.randint(0,9)), str(random.randint(0,9)), str(random.randint(0,9))
    phageID = f'{city}.{one}{two}{three}{four}'
    phageIDs.append(phageID)

prophages_df['phageID'] = phageIDs
prophages_df['phageID'] = prophages_df.apply(lambda row: ''.join([str(row['n']) + '_' + row['phageID']]), axis=1)

prophages_df.to_csv('prophages.txt', sep='\t', index=False)

with open('prophageIDs_used.txt', 'w+') as f:
    f.write('\n'.join(phageIDs))
