%===================文件列表===================%
diagnosis.f90：                 涛师兄cylinder_test3里面的diagnosis，不要改动
diagnosis_gan.f90：             大师姐的diagnosis.f90，不要改动
diagnosis_gan_nc.f90：          大师姐的diagnosis_nc.f90，不要改动
diagnosis_gan_nc_withcom.f90：  大师姐的diagnosis_nc.f90中添加了一些中文注释，不要改动
diagnosis_yang.f90：            涛版diagnosis.f90基础上，将diagnosis_particle改成自己的diagnosis_particle
                                具体改动方案参考大师姐用的diagnosis_nc.f90

%===================使用说明================%
1. 用diagnosis_yang.f90替换原代码中的diagnosis.f90
2. loadparticle.m：读取particle.dat的函数
    output = loadparticle('particle00000000001.dat');
    读出来的粒子数据存储在output.ion和output.ele里面
    output.ion.r：位置，单位m
    output.ion.v：速度，单位m/s
    电子将output.ion换成output.ele即可

%===================可能的问题==================%
1. 输出格式行号可能在程序其他位置重复
    比如有两个1101行
    不过汪老师原代码也有这个问题，应该不影响输出
2. diagnosis_yang和loadparticle.m都没有debug!!!!!!!!!!!!!

