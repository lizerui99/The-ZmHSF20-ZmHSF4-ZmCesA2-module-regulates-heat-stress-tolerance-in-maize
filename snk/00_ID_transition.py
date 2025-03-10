import sys

merge_gtf = sys.argv[1]
input_tab = sys.argv[2]
output_tab = sys.argv[3]

id_dict = dict()
out_file = open(output_tab,'w')

def read_gtf(file_name):
    global id_dict
    with open(file_name,'r') as f:
        for iline in f:
            line = iline.strip()
            if line[0] == "#":
                continue
            else:
                table = line.split('\t')
                if table[2] == "transcript":
                    # gene_id "MSTRG.1"; transcript_id "MSTRG.1.1";
                    # gene_id "MSTRG.15663"; transcript_id "transcript:Zm00001d044095_T001"; ref_gene_id "gene:Zm00001d044095";
                    transcript_table = str(table[-1].replace('"','').replace(' ','')).split(';')
                    #print(transcript_table)
                    gene_id = transcript_table[0][7:]
                    # ['gene_idENSMUSG00000095993', 'transcript_idENSMUST00000179346', 'gene_nameGm21060', 'ref_gene_idENSMUSG00000095993', '']
                    if gene_id in id_dict:
                        if id_dict[gene_id][:5] == 'MSTRG':
                            if transcript_table[-2][:3] == 'ref':
                                ref_id = transcript_table[-2].replace('"', '').replace(' ', '').replace('gene:', '')[
                                         11:]
                            else:
                                ref_id = gene_id
                            id_dict[gene_id] = ref_id
                        else:
                            continue
                    else:
                        if transcript_table[-2][:3] == 'ref':
                            ref_id = transcript_table[-2].replace('"', '').replace(' ', '').replace('gene:', '')[11:]
                        else:
                            ref_id = gene_id
                        id_dict[gene_id] = ref_id



def transition_id(file_name):
    with open(file_name,'r') as f:
        # MSTRG.15663	8.056945333333333	6.866495666666668	9.189581666666667	6.000939	6.999192333333333	6.467595666666667	6.974720666666667	5.406225666666667	6.958373333333333	8.964055	9.860277333333334	4.200996666666667	5.295281666666667	5.5895806666666665	4.241751666666667	5.703395333333333	6.776666666666667	10.229520333333333	3.5721216666666664	8.379552333333335	6.093244333333334	7.593436666666666	8.714449	7.368594333333333	3.2835440000000005	5.805253666666665	9.057698	8.637052333333335	7.705202333333332	7.665607666666666	8.855042333333333	4.800379333333333	4.599570333333333	5.803622999999999	6.536201666666666	6.856888333333334	5.109770999999999	3.7755693333333333	5.399219666666667	5.044272666666667	6.2015	4.5740973333333335	5.289243	4.51007	5.106352333333334	5.5802223333333325	4.830040666666666	5.3569596666666675	9.329452999999999	10.005273	13.161436666666669	6.849592666666666
        # MSTRG.15663	3.886108	9.064914	11.219814	3.571678	8.4792	8.548609	11.170802	11.559008	4.838935	10.268964	4.555315	3.178538	3.798076	8.268828	8.930673	6.020717	6.339352	7.042718	6.750547	5.963434	8.210181	3.209581	3.097803	9.911293	7.257025	3.294721	10.323374	9.175279	9.033858	8.683028	14.824616	10.366987	4.389229	3.55794	3.926698	5.118352	5.327599	5.472244	5.086002	6.057681	7.117104	3.593957	4.17191	5.834146	2.719199	5.998545	7.301748	3.809893	3.740616	9.161426	7.427958	10.055045	10.984987	9.648529	2.01518	2.531677	6.169508	5.062715	10.456805	9.619137	4.113802	8.293477	5.872454	10.245131	7.141787	5.393392	6.299402	14.225379	5.618566	11.948001	4.827503	5.330279	3.188054	3.886208	2.77637	5.872538	5.292724	6.250499	4.457766	11.894598	10.82073	9.470936	4.535185	11.905036	7.966364	4.831453	10.31779	8.933118	5.544599	8.519106	8.288775	8.580903	9.695449	4.318181	5.801237	4.28172	2.666258	5.395142	5.737311	9.316977	4.247619	3.846273	4.04019	4.40964	11.158775	10.073791	6.745657	3.751217	5.237286	4.48757	5.604457	6.01368	2.745604	2.567424	4.189687	3.002651	9.005321	7.518699	1.717696	5.896423	5.760162	7.044909	5.799429	3.603668	8.202663	1.915961	6.968232	2.57303	6.326467	7.867731	2.982177	2.680302	6.096742	4.864375	4.35794	3.486086	9.5941	3.660481	5.164831	3.058201	6.26709	4.19711	4.036733	7.837036	15.533547	8.053785	4.401027	10.28182	8.749889	10.98411	9.763837	18.544457	11.176016	9.022081	2.812602	8.714095
        global id_dict
        num = 1
        for iline in f:
            if num == 1:
                num += 1
                out_file.write(f"{iline}")
            else:
                line = iline.strip().split(',')
                if line[0] in id_dict:
                    line[0] = id_dict[line[0]]
                    str = '\t'.join(line)
                    out_file.write(f"{str}\n")
                else:
                    out_file.write(f"{iline}\n")

read_gtf(merge_gtf)
transition_id(input_tab)

out_file.close()



