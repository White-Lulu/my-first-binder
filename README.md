# Maxwell–Boltzmann Distribution Rates Verified with Two-Dimensional Particle Collisions

## Functionality
- 选择**粒子**(强调“绘图”，而不一定是真的基本粒子)**数目**，**粒子相对分子质量**，**温度**，以及绘图选项（**粒子半径**，**直方图数目**）进行模拟
- 获得粒子速度的直方图分布和拟合曲线，与理论麦克斯韦分布曲线比较

## Tips

- 边界长度和粒子速度的单位都是m(m/s)
- 左图的**粒子颜色**：粒子速度快慢（颜色柱属于左图）；中间图**直方图颜色**：该速度下的粒子数所占总粒子数的比例
- 中间图**红色曲线**：在设定的温度T与相对分子质量M下，理论上的麦克斯韦速率分布曲线；中间图**蓝色曲线**：通过粒子模拟的数据拟合与理论曲线一致的函数形式所得到的（实时）拟合曲线
- 右图**红色折线**：理论曲线值与拟合曲线值的方差；右图**蓝色折线**：拟合曲线值与粒子实际数据值的方差
- 由于动画帧率有限等原因，速度较大的粒子在碰撞时会出现穿模的情况..
- 如果粒子数过少或者粒子半径过大，则数据会与理论有较大出入，导致曲线无法合理拟合，右图两个方差均会出现异常大的值
- 直方图箱数越多，拟合应该越准确

## Ways to run

- **无Jupyter Notebook环境→云运行：https://mybinder.org/v2/gh/White-Lulu/my-first-binder/HEAD?labpath=maxwell-ratesss.ipynb （进入Binder页面后请等待直到进入JupyterLab页面）（可能会加载一会） 对应文件：maxwell-ratesss.ipynb**
  
  ![fdae67628e4397dce5b85ee64fe6af4f](https://github.com/White-Lulu/my-first-binder/assets/173527558/809f1282-e8d9-41c3-a7d2-384da49c643e)

- **在本地Jupyter Notebook运行：maxwell-ratess.ipynb（√多次绘图时不会覆盖上次结果）**

  ![125d6fea52d5227c016277f6627bdf27](https://github.com/White-Lulu/my-first-binder/assets/173527558/c2283c6b-6573-497b-b485-9d8233fd9a8a)

- **在本地Python环境运行：maxwell-rates.py （×没有交互控件，需要手动调参数↓）**
  
  ![image](https://github.com/White-Lulu/my-first-binder/assets/173527558/defc67e0-993c-4cf3-b419-ffa5243ed947)

