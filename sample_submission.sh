#!/bin/sh
#SBATCH -p RM
#SBATCH -t 5:00:00

./eclippe -i /jet/home/sdai/clip_proj/input_data -o /jet/home/sdai/clip_proj/output_dir -t /jet/home/sdai/clip_proj/temp_dir