## Language
- [English](#english)
- [中文](#中文)
---
# English

# Twinkle
A GPU-based High-efficiency binary microlensing code

If you use Twinkle (or part of Twinkle) in your research work, we request citing our paper: Wang, S. et al., 2024, in prep.

If you incorporate Twinkle (or part of Twinkle) in your code, we request specifying it in the README file (or equivalent), including the GitHub link ( https://github.com/AsterLight0626/Twinkle ) and asking the user to cite our paper: Wang, S. et al., 2024, in prep.

# Tutorial
It is recommended to use CUDA for computations. Compile using the provided MakeFile.

In addition to the main function, there are a few parameters specified in `./src/calculation/MacroVariables.h`. These include:

BATCH_SIZE: The number of points sampled on the source edge per iteration
RELTOL: The relative tolerance; computation is completed when (estimated error / magnification) is less than RELTOL
NPOINTS: The maximum number of points that can be sampled on the source edge
NSRCS: The total number of sources the program needs to calculate

## 1. Calculation of Magnification

Use `/src/demo.cpp` as a demonstration; provide relevant parameters to output magnification results.

First, initialize and allocate memory space:

``` 
Twinkle<double> LittleStar;
LittleStar.malloc(0);       // set GPU device number 0
```

Input parameters and time
```
LittleStar.set_path(time,&params);
```

Where, the trajectory parameter $t_0$ is the time when the trajectory is closest to the origin point (mass center), $t_E$ is the time required to traverse the Einstein Radius $\theta_E$, $u_0$ is the closest distance of the trajectory from the origin, $\alpha$ is the angle from the positive x-axis to the trajectory. $\rho$ is the radius of the background source (normalized by $\theta_E$), $q$ is the mass ratio of the two lenses ($0<q<1$, with the larger mass lens on the left), $s$ is the distance between the two lens bodies (normalized by $\theta_E$).

The "time" array with a length of NSRCS, stores the time for each source (corresponding to the observation times). NSRCS is specified in ./src/calculation/MacroVariables.h.


Solve and transfer results from GPU to CPU:
```
LittleStar.solve();
LittleStar.cp_back();
```

Record data:
```
LittleStar.writeto("./");
```

The content recorded is specified in `./src/calculation/init.cpp`. Visualization can be done using programs within the JupyterTools folder.

## 2. Fitting Observed Light Curves
Use `/src/MCMCdemo.cpp` as a demonstration; the MCMC algorithm here is only an example, and it is recommended to replace it with other more mature programs.

The calculation method is consistent with the previous section "Calculation of Magnification". After the calculation is completed, transfer the magnification and error back to the CPU:
```
LittleStar.cp_back_ME();
```

The magnification on the CPU side (host side) is located in LittleStar.h_Mag, whose type is array_t<double>, defined in `./src/calculation/init.h`
Alternatively, use LittleStar.h_Mag.data, which will return the starting pointer of the magnification array.
Similarly, the error data is located in LittleStar.h_Err.

After transferring the magnification and error to the CPU, the MCMC step can be carried out on the CPU side.

## 3. Calculation of Critical Curves and Caustics

Use `/src/caus_demo.cpp` as a demonstration, outputting causcrit.txt, which can be viewed in JupyterTools.

## 4. Use the provided Jupyter notebook for result visualization

Use `./JupyterTools/Visualization.ipynb`

Twinkle's output files are of four types, including critcaus.txt, src_ext.txt, img_pt.txt, and src_pt.txt.

Among them, critcaus.txt stores information about critical curves and caustics, src_ext.txt records simplified information of all calculated points, including magnification and error, etc.

img_pt.txt and src_pt.txt are usually used together, recording the solution details of individual extended sources.

For specific information, see `./src/calculation/init.cpp`

# CPU Version

If you do not have suitable GPU equipment, we also provide [Twinkle_CPU](https://github.com/AsterLight0626/Twinkle_CPU). Its algorithm implementation is identical to Twinkle and will yield exactly the same results.


---
# 中文
# Twinkle
一套基于GPU的高性能二体微引力透镜求解程序

如果您在研究工作中使用了Twinkle（或Twinkle的一部分），请引用我们的文章：Wang, S. et al., 2024, in prep.

如果您将Twinkle（或Twinkle的一部分）整合到您的代码中，请您在README文件（或等效文件）中指明这一点。内容包括GitHub链接 (https://github.com/AsterLight0626/Twinkle)，并要求用户引用我们的论文：Wang, S. et al., 2024, in prep.


# 使用教程 
推荐使用CUDA进行计算。使用提供的 MakeFile 进行编译

在 `main` 函数之外，还有小部分参数在 `./src/calculation/MacroVaribles.h` 中指定。包括：

BATCH_SIZE： 每次迭代在源边缘上采样多少个点

RELTOL：     相对容差，当（估计误差 / 放大率）小于 RELTOL 时完成计算

NPOINTS：    源边缘最多可以采样多少个点

NSRCS：      程序一共需要计算多少个源


## 1.批量计算放大率
使用 `/src/demo.cpp` 做演示，提供相关参数即可输出放大率结果

首先进行初始化，分配内存空间：
``` 
Twinkle<double> LittleStar;
LittleStar.malloc(0);       // set GPU device number 0
```

输入参数和时间：
```
LittleStar.set_path(time,&params);
```
其中，轨迹参数 $t_0$ 是轨迹距离原点（透镜质心）最近的时刻，$t_E$ 是走过 Einstein Radius $\theta_E$ 所需的时间，$u_0$ 是轨迹距离原点的最近距离，$\alpha$ 是轨迹与x轴正方向的夹角。

$\rho$ 是背景源的半径（用 $\theta_E$ 归一化），$q$ 是两个透镜的质量比（$0<q<1$，大质量透镜在左侧），$s$ 是两个透镜体的距离（用 $\theta_E$ 归一化）。

长度为 NSRCS 的 time 列表里存储着每个源的所处时刻（对应观测时刻）。NSRCS 在 `./src/calculation/MacroVaribles.h` 中指定。

求解，并将结果从 GPU 传回 CPU：
```
LittleStar.solve();
LittleStar.cp_back();
```

记录数据：
```
LittleStar.writeto("./");
```
记录的内容在 `./src/calculation/init.cpp` 中指定。可以用 JupyterTools 文件夹内的程序进行可视化。

## 2.拟合观测的光变曲线
使用 `/src/MCMCdemo.cpp` 做演示，此处的 MCMC 算法仅作为示例，建议替换为其他更成熟的程序。

计算方法与上一节“批量计算放大率”一致，计算完成后将放大率和误差传输回CPU：

```
LittleStar.cp_back_ME();
```
在 CPU 端（host 端）的放大率位于 LittleStar.h_Mag。它的类型是 array_t<double>，定义在 `./src/calculation/init.h`
或者使用 LittleStar.h_Mag.data，这将返回放大率数组的开头指针。
类似地，误差的数据位于 LittleStar.h_Err

将放大率和误差传递到CPU后，MCMC步骤在CPU侧进行即可

## 3.计算临界曲线和焦散线
使用 `/src/caus_demo.cpp` 做演示，输出 causcrit.txt，可以在 JupyterTools 中查看。

## 4.使用提供的 jupyter notebook 文件进行结果可视化
使用 `./JupyterTools/Visuallization.ipynb`

Twinkle 的输出文件有四类，包括 critcaus.txt, src_ext.txt, img_pt.txt 和 src_pt.txt

其中 critcaus.txt 存储临界曲线和焦散线的信息，src_ext.txt 记录计算的所有点的简化版信息，包括放大率和误差等。
img_pt.txt 和 src_pt.txt 通常一起使用，记录单个展源的求解细节。
具体信息详见 `./src/calculation/init.cpp`

# CPU 版本
如果你并没有合适的 GPU 设备，我们也提供了 [Twinkle_CPU](https://github.com/AsterLight0626/Twinkle_CPU)。它的算法实现与 Twinkle 完全相同，会给出完全一致的结果。
