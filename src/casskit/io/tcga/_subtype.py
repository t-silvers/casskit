import subprocess

with open("test.tsv", "w") as f:
    for record in SeqIO.parse("/home/fil/Desktop/420_2_03_074.fastq", "fastq"):
        f.write("%s %s %s\n" % (record.id,record.seq, record.format("qual")))


encoding = 'ascii'    # specify the encoding of the CSV data
p2 = subprocess.Popen(["Rscript", "/home/users/tsilvers/casskit/src/casskit/io/tcga/_tcgabiolinks.R"], stdout=subprocess.PIPE)
output = p2.communicate()[0].decode(encoding)
edits = csv.reader(output.splitlines(), delimiter=",")
for row in edits:
    print(row)

p1 = subprocess.Popen(["Rscript", "/home/users/tsilvers/casskit/src/casskit/io/tcga/_tcgabiolinks.R"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)


process = subprocess.run(["Rscript /home/users/tsilvers/casskit/src/casskit/io/tcga/_tcgabiolinks.R"], stdout=subprocess.PIPE, shell=True)
(output, err) = process.communicate()
exit_code = process.wait()


import csv, pprint, subprocess, io

pipe = subprocess.Popen(["Rscript", "/home/users/tsilvers/casskit/src/casskit/io/tcga/_tcgabiolinks.R"], stdout=subprocess.PIPE)
pipeWrapper = io.TextIOWrapper(pipe.stdout)
pipeReader = csv.DictReader(pipeWrapper)
listOfDicts = [ dict(row) for row in pipeReader ]
