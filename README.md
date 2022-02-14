# ESPRH
区域地震目录自动构建流程脚本，整合了开源的各个模块，用于从连续波形自动构建区域地震目录。

Automatic regional-scale earthquake catalog building workflow: EQTransformer + Siamese EQTransformer + PickNet + REAL + HypoInverse.

## 安装说明 Installation
```Bash
conda create -n ESPRH
conda activate ESPRH
conda install python=3.6 tensorflow-gpu=1.14 keras-gpu=2.3.1 h5py=2.10 matplotlib=3.2 pyyaml cudatoolkit cudnn pandas tqdm pyproj jupyter notebook basemap
conda install -c conda-forge obspy
pip install keras-rectified-adam
```
注：只是推理的话不需要GPU也可执行，则将对应的tensorflow-gpu, keras-gpu换成tensorflow, keras。并且不要安装 cudatoolkit和cudnn。

复制该程序包到你的计算机。
Clone this project to your machine. 

```bash
git clone https://github.com/MrXiaoXiao/ESPRH
cd ESPRH
```

REAL和HypoInverse的安装请参照它们对应的说明.HypoInverse在bin下的文件名请设置为‘hyp1.40'.


## 使用说明 Usage
进入目录，依次执行脚本00-06，在default_pipline_config.yaml修改对应参数。

Enter the directory. Execute scripts 00-06. Customize your configurations in file default_pipline_config.yaml.

## 相关的工作 Related research
[1] Wu, Xueshan; Huang, Song; Xiao, Zhuowei; Wang, Yuan (2022): Building Precise Local Submarine Earthquake Catalogs via a Deep-Learning-Empowered Workflow and its Application to the Challenger Deep. Frontiers. Collection. https://doi.org/10.3389/feart.2022.817551 

[2] Shun Yang, Zhuowei Xiao, Yue Zhu, Yumei He, Mingming Jiang, Chit Thet Mon, Guangbing Hou, Myo Thant, Kyaing Sein. (2021, In Preparation). A deep-learning-empowered pipeline for building regional earthquake catalogs and its application to the central Myanmar region.

## 引用 Citation
如果你使用该脚本，请在文章中引用以下工作：

For EQTransformer, please cite:

S. Mostafa Mousavi, William L Ellsworth, Weiqiang Zhu, Lindsay Y Chuang, and Gregory C Beroza. (2020). Earthquake transformer—an attentive deep-learning model for simultaneous earthquake detection and phase picking. Nature Communications 11, 3952. https://doi.org/10.1038/s41467-020-17591-w

For Siamese EQTransformer, please cite:

Zhuowei Xiao, Jian Wang*, Chang Liu, Juan Li, Liang Zhao, and Zhenxing Yao. (2021). Siamese Earthquake Transformer: A pair-input deep-learning model for earthquake detection and phase picking on a seismic array. Journal of Geophysics Research: Solid Earth. https://doi.org/10.1029/2020JB021444

PickNet for phase refinement:

Wang, J., Xiao, Z., Liu, C., Zhao, D., & Yao, Z. (2019). Deep Learning for Picking Seismic Arrival Times. Journal of Geophysical Research: Solid Earth, 124(7), 6612–6624. https://doi.org/10.1029/2019JB017536

REAL for linking seismic phases:

Miao Zhang, William L Ellsworth, and Gregory C Beroza. (2019). Rapid Earthquake Association and Location. Seismological Research Letters, 90(6), 2276–2284. https://doi.org/10.1785/0220190052

HypoInverse for locating earthquakes:

Fred W Klein. (2002). Userʼs Guide to HYPOINVERSE-2000, a Fortran Program to Solve for Earthquake Locations and Magnitudes 4/2002 version. USGS, Open File Report 02-171 Version, 1, 123.

## 近期更新计划
1. 增加02_run_S-EqT步骤的并行加速
2. 增加计算里氏震级（根据震相计算对应于窗内振幅再输入到REAL里面）
(最近事情有点多，大概预计12.20之前更新上）

## 问题反馈 Bug report
如果遇到程序上的问题，请在这个repo开启一个issue（尽量不要邮件联系）。

If you occur any bugs or questions, you can open a new issue in this repo. 

## 邮箱 E-mail
xiaozhuowei@mail.iggcas.ac.cn

## 致谢 Acknowledgments
We would like to thank S. Mostafa Mousavi and his colleagues for developing the EqT model (https://github.com/smousavi05/EQTransformer), which is the base of our S-EqT model.

We would like to thank Miao Zhang for developing REAL (https://github.com/Dal-mzhang/REAL).

We would like to thank Fred Klein for developing HypoInverse (https://www.usgs.gov/software/hypoinverse-earthquake-location)

We would like to thank Yijian Zhou for developing the python interface for HypoInverse (https://github.com/YijianZhou/Hypo-Interface-Py)
