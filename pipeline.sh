#load package
module load cutadapt/2.10
module load STAR/2.7.6a
module load samtools/1.13.0
module load anaconda3/2020.11

#setup anaconda
mv $HOME/.conda $PROJECT/.conda
ln -s $PROJECT/.conda $HOME/.conda

#download package from github
##clipper
git clone https://github.com/YeoLab/clipper.git
cd clipper
conda env create -f environment3.yml
conda activate clipper3
pip install .
##umi_tools
git clone https://github.com/CGATOxford/UMI-tools.git
cd cd UMI-tools
python setup.py install
##fastq_tools
git clone https://github.com/dcjones/fastq-tools.git
bash autogen.sh
./configure
make install
ln -s /usr/bin /jet/home/sdai/bin


#trim adapters