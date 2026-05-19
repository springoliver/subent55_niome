The next phase of the CFTR multi-variant challenge is now LIVE. 🎉 ✨ 
This task extends beyond variant recovery and introduces a clinically relevant objective: predicting therapeutic response from inferred CFTR genotypes.

Challenge Overview
Miners will work with realistic simulated sequencing data representing individuals with:
confirmed Cystic Fibrosis (CF), or
suspected CF-related phenotypes.
The dataset is designed to reflect the complexity of real-world CFTR analysis, including difficult variant combinations and clinically actionable mutations.

What You Will Receive
Simulated diplotype sequencing reads
Reads generated against the GRCh38 human reference assembly
Datasets containing 5000+ CFTR-related variants
Diverse mutation classes, including: SNVs, indels, complex alleles, multi-variant haplotypes
These simulated samples are designed to emulate realistic sequencing and interpretation scenarios encountered in clinical genomics.

Your Objectives
Variant Recovery
Use variant-calling pipelines to identify the variants embedded in the sequencing reads.
This builds directly on the previous task:
alignment
variant calling
genotype reconstruction

Drug Response Prediction (New)
Using the recovered variants and inferred genotypes, predict:
expected CFTR modulator response
likely therapeutic sensitivity
This introduces a translational layer connecting genomic inference to precision medicine.

Validation & Scoring
To ensure fair and biologically meaningful evaluation:
Validation is performed only on variants truly present in the simulated reads
Scoring is weighted by variant difficulty and complexity
Hard-to-detect or clinically significant variants contribute more strongly to the final score
This rewards both:
analytical accuracy
robustness on challenging CFTR architectures

Why This Matters
CFTR is one of the most complex clinically actionable genes due to:
extensive allelic heterogeneity
compound heterozygosity
variant-specific therapeutic effects
This challenge is designed to benchmark methods capable of moving from:
raw sequencing reads → genotype inference → clinical interpretation

Expected Deliverables
Participants should submit:
called variants / inferred genotypes
confidence metrics
predicted drug-response classifications
methodological documentation

Key Focus Areas
Accurate CFTR variant detection
Multi-variant haplotype reconstruction
Clinical interpretation pipelines
Precision medicine inference
Robustness on difficult genomic regions

🔥 Currently, the burn rate is 90%
We will release the remaining emissions after the testing period.

🔢 Updated Scoring Distribution
The new emission distribution for the top 10 miners is now:
[0.2, 0.2, 0.2, 0.1, 0.1, 0.05, 0.05, 0.05, 0.025, 0.025]

⏰ Validator Testing Schedule
During the testing period, validators will run validation every 4 hours at:
00:00, 04:00, 08:00, 12:00, 16:00, 20:00 (UTC)

📈  Leaderboard & Monitoring
The leaderboard is currently updating.
While it's refreshing, you can still monitor validator activity and miner response results through WANDB: 
https://wandb.ai/genomes/niome/table?nw=nwuserjgenome

We'll continue doing our best to contribute to the genomics ecosystem, strengthen the task quality, and support all miners with a transparent and reliable workflow. 🚀 
Thanks for your continued participation and effort across the subnet. 🤝 
