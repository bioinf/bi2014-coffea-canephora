# Read the input data.
# Формат файла "номер гена номер кластера"
clusters = {}
with open('clusters.dmp', 'r') as clusters_data:
    for cluster in clusters_data:
        cluster_num = int(cluster.strip('\n').split(" ")[1])
        gene_num = int(cluster.strip('\n').split(" ")[0])
 
        if cluster_num in clusters.keys():
            clusters[cluster_num].append(gene_num)
        else:
            clusters[cluster_num] = [gene_num]
# На выходе кластер:[гены]
 
# Формат файла "номер гена описание гена"
headers = {}
with open('headers.dmp', 'r') as headers_data:
    headers = { int(header.strip('\n').split(" ")[0]): header.strip('\n')[len(header.strip('\n').split(" ")[0]):]
                for header in headers_data}
 
# Получаем организмы
# {номер гена : организм}
org_names = {}
# организмы
names = set()
for header in headers.items():
    if header[0] in org_names.keys():
        org_names[header[0]].append(header[1][0:4])
    else:
        org_names[header[0]] = [header[1][0:4]]
 
    names.add(header[1][0:4])
 
print(org_names)
 
names = "".join(names).split(' ')[1:]
print (names)
# cluster[1] - номер гена
# cluster[0] - номер кластера к которому он относится
with open('cafe_input.out', 'w') as result_data:
    firstline = "Description" + '\t' + "ID" + '\t'
 
    res_str = ""
    for name in names:
        res_str += name + '_'
        firstline += name + "\t"
 
    firstline = firstline[:-1] + "\n"
    res_str = res_str[:-1]
 
    result_data.writelines(firstline)
 
    print(res_str)
    curr_name_count = {}
 
    for cluster in clusters.items():
        curr_name_count[names[0]] = 0
        curr_name_count[names[1]] = 0
        curr_name_count[names[2]] = 0
        result_line = ""
 
        for gene in cluster[1]:
            if headers[gene].find(names[0]) != -1:
                curr_name_count[names[0]] += 1
            if headers[gene].find(names[1]) != -1:
                 curr_name_count[names[1]] += 1
            if headers[gene].find(names[2]) != -1:
                 curr_name_count[names[2]] += 1
 
        result_line = str("Organisms_order") + '_' + res_str + '\t' + str(cluster[0]) \
                      + '\t' + str(curr_name_count[names[0]]) + '\t' + \
                      str(curr_name_count[names[1]]) + '\t' + str(curr_name_count[names[2]]) + '\n'
 
        result_data.writelines(result_line)