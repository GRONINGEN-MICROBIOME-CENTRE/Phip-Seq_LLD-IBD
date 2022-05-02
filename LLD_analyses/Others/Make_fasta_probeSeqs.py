

Out = "Data/Probes.fa"
with open(Out, "w") as O: pass
with open("Data/df_info_AT_v2_2_SZ_New.csv") as F:
	for line in F:
		l = line.rstrip().split(",")
		if l[0] == "order": continue
		if "(" in l[3]:
			l[3] = l[3].split("(")[0].strip()
		Fasta = ">" + l[0] + ":"+l[4] + "\n" + l[3] + "\n"
		with open(Out, "a") as O:
			O.write(Fasta)
