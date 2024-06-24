import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.optimize import curve_fit

### 参数设置

## 默认

L=2000 # 正方形边界长度 理论单位：m

dt=0.1 # 时间步长 #太大会导致在有限的帧率内穿模（越红越穿模，越紫越正常。。），太小导致视觉速度很慢

s=1000 # 计算步数

alpha_=0.9 # 直方图透明度

k=1.380649*10**(-23) # 玻尔兹曼常数

interval_=100 # 帧之间间隔毫秒数

## 自定义

r=30 # 粒子绘图半径 建议值：L/100

num=300 # 粒子数量

m=29/6.03*10**(-26) # 分子质量 # 数量级和k一样

T=200 # 温度

v_ij=(3*T*k/m)**0.5/1.085 # 初始平均速率=初始均方根速率/1.085  空气分子（m=29）在300k下算出为：508 (m/s) /1.085 = 468 (m/s)

bins_num=100 # 直方图箱数 越多误差越小


### 初始化数据

## 位置数据[均匀分布]
map_x=np.random.rand(num,2)*L

## 速度数据[正态分布]
map_v=np.random.randn(num, 2)*L/10 # 正态分布速度分量
speed=np.linalg.norm(map_v, axis=1) # 速率
mean_speed=np.mean(speed) # 平均速率
scaling_factor=v_ij/mean_speed # 比例因子
map_v*=scaling_factor # 调整速度使平均速率为v_ij
speed=np.linalg.norm(map_v, axis=1) # 速率
mean_speed=np.mean(speed) # 平均速率


### 初始化绘图

## 全局

# 全局图像大小
fig=plt.figure(figsize=(21.8,4.3))

# gridspec分割布局
gs = fig.add_gridspec(1, 5, width_ratios=[7,1.5,6,0.3,6])

# 定义绘图区域到分割的区域
ax1 = fig.add_subplot(gs[0])
ax2 = fig.add_subplot(gs[2])
ax3 = fig.add_subplot(gs[4])

## 粒子运动图

# 坐标轴等于边界长
ax1.set_xlim(0,L) # 理论单位：m
ax1.set_ylim(0,L)

# 散点图，粒子大小，速率与颜色映射，颜色主题
scat=ax1.scatter(map_x[:,0],map_x[:,1],s=(r*25*10/L)**2*3.1416,c=speed,cmap='rainbow')

# 颜色柱
norm=plt.Normalize(vmin=np.min(speed), vmax=np.max(speed))
cbar=plt.colorbar(scat,ax=ax1,fraction=0.046,pad=0.04,norm=norm)
cbar.set_label('Speed')

# 标题
ax1.set_title('Particle Collision Simulation')

## 粒子数据直方图

# bins：箱数 density：是否归一化 alpha：透明度
ax2.hist(speed,bins=bins_num,density=True,label='Experimental')

# 网格线
ax2.grid(alpha=0.15)

## 理论函数曲线

v=np.linspace(0,np.max(speed),100) #自变量坐标 计算点数目
dn_N=2**0.5*3.1416**(-0.5)*3**1.5*v_ij**(-3)*np.exp(-v**2*1.5/v_ij**2)*(v**2) #这里应该用初始速度计算

# 绘图，注释，图例
ax2.plot(v, dn_N, 'r', lw=2,label='Analytical') #lw：线条宽度

## 拟合函数定义

def maxwell_boltzmann(v, a): #curvefit()需要函数输入参数为(自变量，参数1，...)
    return (4/np.sqrt(np.pi))*(a**1.5)*np.exp(-a*v**2)*v**2

## 方差图

# 存储loss的列表
data_fitting_loss=[]
fitting_analytical_loss=[]
frames=[]

# 绘图，网格线，注释，图例
loss1_data,=ax3.plot([],[],'b-',label='Data vs Fitting') #逗号解包语法从ax3.plot返回的列表中提取出Line2D对象
loss2_data,=ax3.plot([],[],'r-',label='Fitting vs Analytical')
ax3.grid(True, color='gray', alpha=0.3)
ax3.set_title('Loss Over Time')
ax3.set_xlabel('Frame')
ax3.set_ylabel('Loss')


###碰撞处理

## 边界碰撞

def wall_collisions(map_x,map_v, L, r):
    for i in range(2):
        mask=map_x[:,i]<r #碰到左/下边界
        map_x[mask,i]=r #退回到位置r
        map_v[mask,i]=-map_v[mask,i] #对称速度反向

        mask=map_x[:,i]>L-r #碰到右/上边界
        map_x[mask,i]=L-r #退回到位置L-r
        map_v[mask,i]=-map_v[mask,i]
    return map_x,map_v

## 粒子间碰撞

def particle_collisions(map_x,map_v,r):
    for i in range(num):
        for j in range(i+1,num):
            xy_ij=map_x[i]-map_x[j] # 数组相减 不重复，计算所有粒子间各距离
            d_ij=np.linalg.norm(xy_ij) # 线性代数，求范数 矩阵整体元素平方和开根号，不保留矩阵二维特性
            if d_ij<=2*r:  # 如果碰撞 # 因为粒子相同，所以交换速度
                xita=xy_ij/d_ij # 碰撞方向关于坐标轴的角度 cos，sin
                vxy_ij=map_v[i]-map_v[j] # 两个粒子的速度分量差
                vxy_ij2=np.dot(vxy_ij,xita) # dot 矩阵乘法 得到速度差分量在碰撞方向的投影（速度只在碰撞方向交换，在垂直碰撞方向不变）
                if vxy_ij2<0: # 通过去掉差值实现速度交换
                    map_v[i]-=vxy_ij2*xita # 速度大的减去差值
                    map_v[j]+=vxy_ij2*xita # 速度小的增加差值
    return map_v


### 计算，更新，绘图

def update(frame):

    global map_x,map_v
    
    ## 计算粒子数据

    map_x+=map_v*dt # 根据速度得到下一步位置
    map_x,map_v=wall_collisions(map_x,map_v,L,r) # 处理边界碰撞
    map_v=particle_collisions(map_x,map_v,r) # 处理粒子间碰撞
    speed=np.linalg.norm(map_v,axis=1) # 算速率

    ## 绘制粒子运动图

    scat.set_offsets(map_x)
    scat.set_array(speed) # 更新粒子颜色
    
    ## ax2由于需要用矩形填色，不能使用offsets，故每次删掉整图，再重画

    ax2.cla()

    ## 计算直方图

    counts,bins,patches=ax2.hist(speed, bins=bins_num,density=True,color='g',alpha=alpha_,label='Experimental')
    ax2.grid(True, zorder=0,color='gray',alpha=0.15)
    bin_colors=plt.cm.rainbow(counts/max(counts))  # 根据占比数目设置颜色 cmap的数字映射为[0,1] max(counts)作为颜色1，则其他颜色映射应按其比例(均匀映射)
    
    ## 绘制直方图

    for count,color,bin_start,bin_end in zip(counts,bin_colors,bins[:-1],bins[1:]): 
        ax2.fill_between([bin_start,bin_end],0,count,color=color) # 对每个bin矩形内填色，通过错位相减确定底边区域
    
    ## 绘制理论函数曲线

    ax2.plot(v,dn_N,'r',lw=2,label='Analytical')
    ax2.set_title('Speed Distribution')
    ax2.set_xlabel('Speed')
    ax2.set_ylabel('Probability Density')

    ## 计算，绘制拟合曲线

    bin_centers=0.5*(bins[1:]+bins[:-1]) # 每个bin的中点
    initial_guess=[0.0001] # 初始参数猜测（出现无法拟合的情况）
    popt=curve_fit(maxwell_boltzmann, bin_centers, counts, p0=initial_guess, bounds=(0, np.inf))[0]# curvefit:非线性最小二乘法拟合 bound：参数边界
    a_fit=popt[0] # curvefit拟合的参数a # pcov: popt的估计协方差
    dn_fit=maxwell_boltzmann(v,a_fit) # 计算拟合函数
    ax2.plot(v,dn_fit,'b',lw=2,label='Fitting') # 绘图

    ## 计算loss

    data_fitting=np.sum((counts-maxwell_boltzmann(bin_centers,a_fit))**2) # 实验与拟合的方差
    fitting_analytical=np.sum((maxwell_boltzmann(v,a_fit)-dn_N)**2) # 拟合与理论的方差

    ## 更新loss

    frames.append(frame)
    data_fitting_loss.append(data_fitting)
    fitting_analytical_loss.append(fitting_analytical)

    ## 绘制loss

    loss1_data.set_data(frames,data_fitting_loss)
    loss2_data.set_data(frames,fitting_analytical_loss)
    ax3.set_xlim(0,s)
    ax3.set_ylim(0,max(max(data_fitting_loss),max(fitting_analytical_loss))*1.1) # 动态改变y轴尺度 适当留一些
    ax3.legend()


### 运行

ani=animation.FuncAnimation(fig,update,frames=s,interval=interval_)
plt.show()
