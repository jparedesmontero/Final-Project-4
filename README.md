
# Investigating Antibiotic Resistance Genes in Healthcare-Associated Infections through Sequence Alignment 
Problem Statement
# Biological Question
Does Pseudomonas aeruginosa contain genes associated with antibiotic resistance, particularly beta-lactamases?
# Hypothesis
The bacterial genome contains genes that show significant similarity to known antibiotic resistance genes, indicating potential resistance to specific antibiotics.
# Significance
Understanding antibiotic resistance mechanisms is essential for public health, as resistant bacteria pose significant challenges to treatment. Identifying resistance genes helps in monitoring the spread of resistance and developing strategies for antibiotic stewardship.
# Goal
Determine whether Pseudomonas aeruginosa harbors genes homologous to known antibiotic resistance genes and assess their potential impact on resistance phenotypes.


SRA numbers - 100 sample size for project
accession list created

# Tools and Software
   #FastQC: Performs quality control checks on raw sequencing data, generating reports about read quality, GC content, adapter content, etc.

#BLAST https://blast.ncbi.nlm.nih.gov/Blast.cgi 
    A widely used tool for sequence alignment. Use the blastn or blastx command to align your reads to a reference database.
    
#ResFinder Database: A specialized database for antibiotic resistance genes.http://genepi.food.dtu.dk/resfinder
    
#CARD https://www.google.com/url?sa=t&source=web&rct=j&opi=89978449&url=https://card.mcmaster.ca/&ved=2ahUKEwiNr5LW-6qMAxWXk4kEHam0HgkQFnoECAkQAQ&usg=AOvVaw3ZCtqK2lNkJ0JlhWV_bZFc

SLURM

SRA-TOOLKIT

****Need to download the SRA files onto the HBC yet. 


### Assemble beta-lactamase from 100 samples from pennsylvania hospitals
- Group successfully assembled the beta-lactamase gene and produced bam files with index bai.
- We will use the bam files to determine whether the sample has the gene or not.
- It appears that the beta-lactamase gene was mapped in all of the samples. To confirm, you can run the followinf script
- Create script:
```
vi presence.sh
```
- Type I to edit file and paste the following script:
```
#!/bin/bash

# Output file
echo "Sample,Gene_Present" > gene_presence.csv

# Loop through each BAM file
for bamfile in SRR*.bam; do
    sample=$(basename "$bamfile" .bam)
    
    # Count mapped reads
    mapped=$(samtools view -c -F 4 "$bamfile")
    
    # Determine gene presence (1 = present, 0 = absent)
    if [ "$mapped" -gt 0 ]; then
        echo "$sample,1" >> gene_presence.csv
    else
        echo "$sample,0" >> gene_presence.csv
    fi
done
```
- Inspct the `gene_presence.csv` file to confirm the beta-lactamase gene was successfully mapped in all samples.


## 1. Call SNPs to study the gene variation 
- Create snpcall.slurm
```
vi snpcall.slurm
```
- Type I to edit and paste (make sure you add your email):
```
#!/bin/bash
#SBATCH --job-name=call_snps
#SBATCH --array=0-99         # Adjust based on number of samples (0 to N-1)
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=04:00:00
#SBATCH --output=call_snps_%A_%a.out
#SBATCH --error=call_snps_%A_%a.err
#SBATCH --mail-user=your.email@svsu.edu
#SBATCH --mail-type=ALL

# Load required modules
module load samtools
module load bcftools

# Index reference (only needed once, but safe here)
if [ ! -f reference.fasta.fai ]; then
    samtools faidx reference.fasta
fi

# Create output directory if not exists
mkdir -p vcfs

# Get all BAM files that start with SRR
BAM_FILES=(SRR*.bam)

# Get the filename for this SLURM array task
BAM_FILE=${BAM_FILES[$SLURM_ARRAY_TASK_ID]}
SAMPLE=$(basename "$BAM_FILE" .bam)

# Call variants, output uncompressed VCF
bcftools mpileup -f reference.fasta "$BAM_FILE" | \
bcftools call -mv -Ov -o "vcfs/${SAMPLE}.vcf"

# Compress and index the VCF
bgzip -c "vcfs/${SAMPLE}.vcf" > "vcfs/${SAMPLE}.vcf.gz"
tabix -p vcf "vcfs/${SAMPLE}.vcf.gz"
```
- Run script
```
sbatch snpcall.slurm
```
- vcf files wil be stored in folder labeled vcfs

# VISUALIZATION
## Mutational Load Per Sample
- Create python script
```
vi mutational_load.py
```
- Type I to edit:
```Python
import os
import pandas as pd
import matplotlib.pyplot as plt

vcf_dir = "vcfs"  # Path to your VCF files
snp_counts = {}

# Loop through each uncompressed VCF file
for file in os.listdir(vcf_dir):
    if file.endswith(".vcf") and not file.endswith(".vcf.gz"):
        sample_id = file.replace(".vcf", "")
        count = 0
        with open(os.path.join(vcf_dir, file), "r") as f:
            for line in f:
                if not line.startswith("#"):
                    count += 1
        snp_counts[sample_id] = count

# Convert to DataFrame
df = pd.DataFrame(list(snp_counts.items()), columns=["Sample", "SNP_Count"])
df = df.sort_values("SNP_Count", ascending=False)

# Plot
plt.figure(figsize=(14, 6))
plt.bar(df["Sample"], df["SNP_Count"], color="royalblue")
plt.axhline(df["SNP_Count"].mean(), color="red", linestyle="--", label="Mean SNPs")
plt.xticks(rotation=90, fontsize=6)
plt.ylabel("Number of SNPs")
plt.xlabel("Sample ID")
plt.title("Mutational Load in Beta-lactamase Gene Across Isolates")
plt.legend()
plt.tight_layout()
plt.savefig("mutational_load_barplot.png", dpi=300)
plt.show()
```
- Run python script
```
module load python/3.8.6
python mutational_load.py
```
- Push and Inspect mutational_load_barplot.png, if good, to poster

### IF PYTHON GIVES YOU TROUBLE YOU CAN PLOT IN R
- Create script
```
vi snps_count.sh
```
- Type I to edit:
```
#!/bin/bash

echo "Sample,SNP_Count" > snp_counts.csv

for vcf in vcfs/*.vcf; do
    sample=$(basename "$vcf" .vcf)
    count=$(grep -vc '^#' "$vcf")
    echo "$sample,$count" >> snp_counts.csv
done
```
- run it
```
bash snps_count.sh
```
- Push csv file to repo and load it to R to build grah
- In Rstudio, the following script can be used
```R
# Read the CSV file
df <- read.csv("snp_counts.csv")

# Basic bar plot
barplot(
  df$SNP_Count,
  names.arg = df$Sample,
  las = 2,                      # Rotate sample labels
  col = "steelblue",
  main = "Mutational Load in Beta-lactamase Gene Across Isolates",
  ylab = "Number of SNPs",
  ylim = c(0, max(df$SNP_Count) + 2)
)

# Add a dashed red line for the mean
abline(h = mean(df$SNP_Count), col = "red", lty = 2)

# Add legend
legend("topright", legend = "Mean SNPs", col = "red", lty = 2, bty = "n")
```

## Mutation Hotspots in the betalactamase gene
- Create python script
```
vi hotspot.py
```
- Type I to edit and paste:
```
import os
import pandas as pd
import matplotlib.pyplot as plt

vcf_dir = "vcfs"  # Update if needed
position_counts = {}

for file in os.listdir(vcf_dir):
    if file.endswith(".vcf") and not file.endswith(".vcf.gz"):
        with open(os.path.join(vcf_dir, file), "r") as f:
            for line in f:
                if not line.startswith("#"):
                    pos = int(line.strip().split("\t")[1])
                    position_counts[pos] = position_counts.get(pos, 0) + 1

df = pd.DataFrame(sorted(position_counts.items()), columns=["Position", "Frequency"])

plt.figure(figsize=(12, 5))
plt.stem(df["Position"], df["Frequency"], basefmt=" ", linefmt='gray', markerfmt='ro')
plt.scatter(df["Position"], df["Frequency"], color='red', s=60)
plt.title("Mutation Hotspots in the Beta-lactamase Gene Suggest Functional Adaptation")
plt.xlabel("Position in Gene")
plt.ylabel("Number of Samples with Mutation")
plt.grid(axis='y', linestyle='--', alpha=0.6)
plt.tight_layout()
plt.savefig("mutation_hotspots_lollipop.png", dpi=300)
plt.show()
```
- Run python script
```
module load python/3.8.6
python hotspot.py
```
- Push and Inspect mutation_hotspots_lollipop.png, if good, to poster

### IF PYTHON GIVES YOU TROUBLE YOU CAN PLOT IN R
- Create script
```
vi hotsopt.sh
```
- Type I to edit:
```
#!/bin/bash

# Temporary file with all sample-position pairs
echo "Sample,Position" > mutation_positions.csv

for vcf in vcfs/*.vcf; do
    sample=$(basename "$vcf" .vcf)
    grep -v '^#' "$vcf" | cut -f2 | awk -v s="$sample" '{print s","$1}' >> mutation_positions.csv
done

# Now summarize how many unique samples have a mutation at each position
# Output: Position,Frequency
tail -n +2 mutation_positions.csv | sort -u | cut -d',' -f2 | sort | uniq -c | awk '{print $2","$1}' > mutation_hotspot_counts.csv

# Add header
sed -i '1iPosition,Frequency' mutation_hotspot_counts.csv

```
- run it
```
bash hotspot.sh
```
- Push mutation_hotspot_counts.csv file to repo and load it to R to build grah
- In Rstudio, the following script can be used
```R
df <- read.csv("mutation_hotspot_counts.csv")

ggplot(df, aes(x = Position, y = Frequency)) +
  geom_segment(aes(xend = Position, yend = 0), color = "gray") +
  geom_point(color = "red", size = 3) +
  labs(
    title = "Mutation Hotspots in the Beta-lactamase Gene Suggest Functional Adaptation",
    x = "Position in Gene",
    y = "Number of Samples with SNP"
  ) +
  theme_minimal()
```

