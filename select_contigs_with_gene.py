
spl = ["Fucus-distichus", "Heribaudiella-fluviatilis", "Himanthalia-elongata", "Macrocystis-pyrifera_MALE",
       "Pelvetia-canaliculata", "Saccorhiza-dermatodea", "Saccorhiza-polyschides_MALE", "Desmarestia-dudresnayi"]

spl2 = ["Desmarestia-dudresnayi"]

for sp in spl:
    gff_path = f"data/Genomes/{sp}/{sp}.gff"
    fa_path = f"data/Genomes/{sp}/{sp}.fa"
    new_fa = f"data/Genomes/{sp}/new_{sp}.fa"
    contigs_list = set()
    with open(gff_path, "r") as gff:
        for l in gff:
            if "ID" in l:
                l = l.split("ID=")
                id = l[1].split(";")[0][5:].split(".")[0]
                contigs_list.add(id)
        print(f"{sp} : {len(contigs_list)}")
    with open(fa_path, "r") as ifa, open(new_fa, "w") as ofa:
        write = False
        for l in ifa:
            if l[0] == ">":
                write = False
                if l[1:-1] in contigs_list:
                    write = True
            if write:
                ofa.write(l)
