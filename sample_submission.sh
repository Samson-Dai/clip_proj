#!/bin/sh
#SBATCH -p RM-shared
#SBATCH -t 5:00:00
#SBATCH --mem=2000M

./eclippe -i /jet/home/sdai/clip_proj/input_data -o /jet/home/sdai/clip_proj/output_dir -t /jet/home/sdai/clip_proj/temp_dir