#load package
module load cutadapt/2.10
module load STAR/2.7.6a
module load samtools/1.13.0
module load anaconda3/2020.11

#setup anaconda
mv $HOME/.conda $PROJECT/.conda
ln -s $PROJECT/.conda $HOME/.conda

#download package via anaconda
conda install -c bioconda -c conda-forge umi_tools
conda install -c bioconda fastq-tools

#download package from github
#https://github.com/YeoLab/clipper
git clone https://github.com/YeoLab/clipper.git
cd clipper
conda env create -f environment3.yml
conda activate clipper3
pip install .


#trim adapters