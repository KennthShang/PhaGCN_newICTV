#!/bin/bash
#SBATCH -o log.out
#SBATCH -J Rscript
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=50G
#SBATCH --time=72:00:00
#SBATCH --export=all
#SBATCH --gres=gpu:1






time python run_Speed_up.py --contigs test.fasta --len 1999
#python tmp.py
#bash code/compress_script.sh
#time python cnn_train.py
#time python cnn_preprocessing.py

#time blastp -query test_neg.fa -db db/db -outfmt 6 -evalue 1e-5 -num_threads 8 -out "test_neg.tab" -max_target_seqs 1
#time blastp -query pvp_v3.fa -db db/db -outfmt 6 -evalue 1e-200 -num_threads 32 -out "pvp.tab"
#time python smote.py

#python preprocessing_kmer.py
#python preprocessing.py --file1 pvpred/test_pos --file2 pvpred/test_neg
#python preprocessing.py --file1 pvpred/train_pos --file2 pvpred/train_neg
#time python preprocessing.py --file1 pvpred/smote_train_pos --file2 pvpred/smote_train_neg

#Rscript generate_cgr_with_aaindex.R train_neg AMINO
#Rscript generate_cgr_with_aaindex.R train_pos AMINO


#python tmp.py
#python evenscale.py
#Rscript generate_cgr_with_fft.R train_pos AMINO
#Rscript generate_cgr_with_fft.R train_neg AMINO
#Rscript generate_cgr_with_fft.R test_neg AMINO
#Rscript generate_cgr_with_fft.R test_pos AMINO
#Rscript generate_cgr_with_aaindex_update.R test_pos AMINO
#Rscript generate_cgr_with_aaindex_update.R test_neg AMINO
#python bert_pred.py
#torchrun --nproc_per_node=4 train.py
#python finetune.py
#python -m torch.distributed.launch --nproc_per_node=4 finetune.py
