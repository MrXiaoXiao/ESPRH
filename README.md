# ESPRH
区域地震目录自动构建流程脚本，整合了开源的各个模块，用于从连续波形自动构建区域地震目录。

Automatic regional-scale earthquake catalog building workflow: EQTransformer + Siamese EQTransformer + PickNet + REAL + HypoInverse.

## 安装说明 Installation
```Bash
conda create -n seqt
conda activate seqt
conda install python=3.6 tensorflow-gpu=1.14 keras-gpu=2.3.1 h5py=2.10 matplotlib=3.2 pyyaml cudatoolkit cudnn pandas tqdm pyproj jupyter notebook basemap
conda install -c conda-forge obspy
pip install keras-rectified-adam
```
注：只是推理的话不需要GPU也可执行，则将对应的tensorflow-gpu, keras-gpu换成tensorflow, keras。并且不要安装 cudatoolkit和cudnn。

## 使用说明 Usage
解压文件然后进入目录，依次执行脚本00-06，在default_pipline_config.yaml修改对应参数。
Unzip the package and enter the directory. Execute scripts 00-06. Customize your configurations in file default_pipline_config.yaml.

## 应用该流程的工作 Related research
[1] Xueshan Wu, Song Huang, Zhuowei Xiao, Yuan Wang. (2021, Under Reveiw). Building precise local submarine earthquake catalogs via a deep-learning-empowered workflow and its application to the Challenger Deep.
[2] Shun Yang, Zhuowei Xiao, Yue Zhu, Yumei He, Mingming Jiang, Chit Thet Mon, Guangbing Hou, Myo Thant, Kyaing Sein. (2021, In Preparation). A deep-learning-empowered pipeline for building regional earthquake catalogs and its application to the central Myanmar region.

## 引用 Citation
如果你使用改脚本，请在文章中引用以下工作：
If you use the S-EqT codes in your research, please cite:

Zhuowei Xiao, Jian Wang*, Chang Liu, Juan Li, Liang Zhao, and Zhenxing Yao. (2021). Siamese Earthquake Transformer: A pair-input deep-learning model for earthquake detection and phase picking on a seismic array. Journal of Geophysics Research: Solid Earth. https://doi.org/10.1029/2020JB021444

and

S. Mostafa Mousavi, William L Ellsworth, Weiqiang Zhu, Lindsay Y Chuang, and Gregory C Beroza. (2020). Earthquake transformer—an attentive deep-learning model for simultaneous earthquake detection and phase picking. Nature Communications 11, 3952. https://doi.org/10.1038/s41467-020-17591-w

PickNet for phase refinement:

Wang, J.*, Xiao, Z.*, Liu, C., Zhao, D., & Yao, Z. (2019). Deep Learning for Picking Seismic Arrival Times. Journal of Geophysical Research: Solid Earth, 124(7), 6612–6624. https://doi.org/10.1029/2019JB017536

REAL for linking seismic phases:

Miao Zhang, William L Ellsworth, and Gregory C Beroza. (2019). Rapid Earthquake Association and Location. Seismological Research Letters, 90(6), 2276–2284. https://doi.org/10.1785/0220190052

HypoInverse for locating earthquakes:

Fred W Klein. (2002). Userʼs Guide to HYPOINVERSE-2000, a Fortran Program to Solve for Earthquake Locations and Magnitudes 4/2002 version. USGS, Open File Report 02-171 Version, 1, 123.

## 问题反馈 Bug report
如果遇到程序上的问题，请在这个repo开启一个issue（尽量不要邮件联系）。
If you occur any bugs or questions, you can open a new issue in this repo. 

## 致谢 Acknowledgments
We would like to thank S. Mostafa Mousavi and his colleagues for developing the EqT model (https://github.com/smousavi05/EQTransformer), which is the base of our S-EqT model.

We would like to thank Miao Zhang for developing REAL (https://github.com/Dal-mzhang/REAL).

We would like to thank Fred Klein for developing HypoInverse (https://www.usgs.gov/software/hypoinverse-earthquake-location)

We would like to thank Yijian Zhou for developing the python interface for HypoInverse (https://github.com/YijianZhou/Hypo-Interface-Py)
