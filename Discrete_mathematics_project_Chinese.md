<h1 style="text-align: center;">对DNA计算解决布尔可满足性问题的综述与代码仿真</h1>
  
<h3 style="text-align: center;">蒋宏骏<sup>1</sup></h3>
  
<p style="text-align: center;">
  <sup>1</sup>南方科技大学，广东深圳 <br/>
  E-mail: 12411912@mail.sustech.edu.cn
</p>
  
***
  
## 一 摘要
  
NP完全问题、布尔可满足性问题（SAT）是逻辑学、计算机科学的核心问题，其解决方案可影响到数类优化问题的解决效率，因此一直是学者们的关注对象。DNA计算是一种主要利用DNA等生物分子的巨大并行性、高密度存储能力和可设计杂交/脱离操作来执行计算的新型范式，为SAT问题提供了良好的实验平台。本文先回顾 NP问题、SAT/3-SAT问题及Cook-Levin归约思路；随后主要介绍 DNA 计算的核心原理，并选取 *Qinghua Liu(2000)*[^1] 和 *Ravinderjit S. Braich等(2002)*[^2] 两项具有代表性的DNA计算解决SAT问题的实例，剖析其原理、细节以及后续技术路线发展；随后结合python仿真代码解析DNA计算中编码、筛选、放大与读出的步骤，利用代码片段演示“扩增-筛选-连接-电泳-读出”等生化抽象的数字化实现，并添加可视化展示DNA演化流程；最后从误差、复杂度、可编程性与未来路线展望 DNA 计算在更大规模SAT问题与其他复杂问题上的可能后续发展路径。
  
**关键词**：DNA计算，SAT问题、并行计算、生化反应仿真、生化模拟器
  
***
  
## 二 引言
  
布尔可满足性问题（SAT）是逻辑学与计算机科学的核心问题，旨在判断一个给定的合取范式布尔公式是否存在一组变量赋值使其值为真。它也是第一个被证明的NP完全问题，存在从任意 <img src="https://latex.codecogs.com/gif.latex?L\in\text{NP}"/> 的多项式归约 <img src="https://latex.codecogs.com/gif.latex?L%20\le_P%20\text{SAT}"/>，这意味着任何NP问题都可以在多项式时间内转化为SAT问题。因此，其求解效率直接关系到许多组合优化问题。3-SAT限制每个子句最多三个文字，仍保持 NP-完全性。经典求解器（DPLL/CDCL、局部搜索、Survey Propagation）在许多实例上高效，但在最坏情况（相变区域）仍面临指数时间：最坏情形下已知上界 <img src="https://latex.codecogs.com/gif.latex?O^*(1.439^n)"/>、下界 <img src="https://latex.codecogs.com/gif.latex?\Omega(1.3^n)"/>。DNA 计算以分子并行为特性，能够一次性覆盖极大规模候选解，通过杂交选择、酶切、PCR 扩增、分离和读出完成生成—筛选—放大—读出流程，为 NP 问题提供了高效的并行路径。*Liu* 在2000年完成的实验是首个在固体表面（金膜）上演示SAT计算的范例，解决了一个4变量、4子句的3-SAT问题。*Ravinderjit S. Braich等*在2002年采用了“贴纸模型”，实验解决了当时规模空前的20变量、24子句的3-SAT问题。虽然DNA计算受到诸多实际限制，如电泳速率慢、提前制备复杂等，但其并行性理念与碱基互补配对的编码范式被诸多计算机算法借鉴，如图像神经网络（GNN）、DNA遗传算法等。本文目标：在数学与实验两个尺度上链接 SAT 与 DNA 计算，并给出可复现的 Python 模拟器示例与可视化DNA演化过程，将公式—生化操作—代码接口三者对应起来。
  
***
  
## 三 背景理论
### 3.1 NP与NP完全
NP问题是指一类决策问题，其解的正确性能在多项式时间内被验证。一个语言 <img src="https://latex.codecogs.com/gif.latex?L\subseteq\{0,1\}^*"/> 属于NP类，当且仅当存在多项式 <img src="https://latex.codecogs.com/gif.latex?p"/> 与多项式时间判定器（图灵机） <img src="https://latex.codecogs.com/gif.latex?V"/>，使得
  
<p align="center"><img src="https://latex.codecogs.com/gif.latex?x\in%20L%20\Leftrightarrow%20\exists%20y,\%20|y|\le%20p(|x|),\%20V(x,y)=1\quad"/></p>  
  
  
问题 <img src="https://latex.codecogs.com/gif.latex?A"/> 为 NP 完全,需满足：属于NP类，且对于任意问题 <img src="https://latex.codecogs.com/gif.latex?B\in\text{NP}"/>，都存在一个多项式时间归约，即：
<p align="center"><img src="https://latex.codecogs.com/gif.latex?(1)%20A\in\text{NP}；(2)%20\forall%20B\in\text{NP},\%20B\le_P%20A。"/></p>  
  
  
NP完全问题是NP类中满足"任何NP问题都可以在多项式时间内归约到该问题"性质的问题。这类问题具有重要的理论意义，如果能找到任何一个NP完全问题的多项式时间算法，也许就可以推导出P=NP这一重要结论。
  
### 3.2 SAT定义与Cook-Levin归约
  
布尔可满足性问题（SAT）定义如下：  
给定一个布尔公式 <img src="https://latex.codecogs.com/gif.latex?\phi(x_1,x_2,\dots,x_n)"/>，其合取范式（CNF）形为
  
<p align="center"><img src="https://latex.codecogs.com/gif.latex?\phi(x_1,\dots,x_n)=\bigwedge_{i=1}^m%20C_i,\quad%20C_i=\bigvee_{j\in%20S_i}%20l_j,\%20l_j\in\{x_k,\neg%20x_k\}"/></p>  
  
  
判定：是否存在赋值 <img src="https://latex.codecogs.com/gif.latex?\sigma:\{x_i\}\to\{\text{T},\text{F}\}"/> 使 <img src="https://latex.codecogs.com/gif.latex?\phi(\sigma)=\text{T}"/>。
Cook-Levin 归约：用布尔变量编码图灵机空间-时间-状态-符号轨迹，构造以下约束：
  
<div style="text-align:center">
  
  |约束|描述|
  |:-- |:-- |
  |唯一性|每时间步恰一状态、每格恰一符号|
  |合法转移|相邻时间步满足转移函数|
  |初始配置|输入与起始状态一致|
  |接受条件|存在时间步处于接受状态|
  
</div>
  
  以此为基础进行进一步构造，可证得 <img src="https://latex.codecogs.com/gif.latex?\phi_x"/> 可满足 <img src="https://latex.codecogs.com/gif.latex?\Leftrightarrow"/> <img src="https://latex.codecogs.com/gif.latex?x"/> 。
  
  
### 3.3 3-SAT问题及其传统求解困难
  
3-SAT 是SAT的限制版本：每个子句恰好包含三个文字,仍是NP完全的 *（1972，Karp）*。传统硅基计算机求解3-SAT的主要算法包括回溯法、贪心算法、局部搜索、启发式随机行走算法。尽管现代的SAT求解器在实例上表现不错，但在最坏情况下仍会退化到指数时间。已知下界为 <img src="https://latex.codecogs.com/gif.latex?\Omega(1.307^n)"/>（2023年），而上界为 <img src="https://latex.codecogs.com/gif.latex?O^*(1.439^n)"/>。对 <img src="https://latex.codecogs.com/gif.latex?n%20\ge%201000"/> 的随机3-SAT实例来说，求解时间仍不可控。这促使人们寻找非冯·诺伊曼架构之外的计算范式，DNA计算正是其中最具代表性的非传统计算模型之一。
  
### 3.4 DNA计算模型及其在SAT问题上的应用能力
DNA计算的建模思想是将计算问题的某些操作或要素映射到DNA分子及碱基片段上。其核心思想是利用DNA碱基配对原则（A-T, G-C）来编码信息，并通过生物化学操作实现计算过程。基础的DNA计算模型包括以下两个步骤：
  
- 编码策略：使用DNA序列编码问题的变量、约束条件和候选解信息，如每个变量对应一个碱基对、子句对应互补序列、候选解对应一个序列片段或DNA单链。
  
- 生化操作：利用一系列生化操作来模拟计算过程中的逻辑运算和筛选，包括：
  - 杂交——通过互补片段结合实现逻辑匹配；
  - 酶切——使用限制性内切/外切酶去除不满足条件的序列；
  - PCR扩增——复制特定DNA片段以增加浓度；
  - 电泳分离——根据长度或质量分离不同的DNA分子。
  
以一个简单的SAT问题为例:假设有布尔变量<img src="https://latex.codecogs.com/gif.latex?x_1,%20x_2,%20x_3"/>，需要找到使 <img src="https://latex.codecogs.com/gif.latex?(x_1%20\vee%20\bar%20x_2)%20\wedge%20(\bar%20x_1%20\vee%20x_3)"/> 为真的赋值。（1）编码阶段：为每个变量创建代表真假的DNA序列，比如，对于变量<img src="https://latex.codecogs.com/gif.latex?x_1"/>，可以用序列ATCG代表<img src="https://latex.codecogs.com/gif.latex?x_1=true"/>，用GCAT代表<img src="https://latex.codecogs.com/gif.latex?x_1=false"/>。（2）并行生成：通过在溶液环境中混合所有可能的组合，同时生成所有<img src="https://latex.codecogs.com/gif.latex?2^3=8"/>种可能的赋值组合，即，使得所有待选解能被同时验证。（3）筛选阶段：针对每个子句设计相应的探针，通过杂交互补片段确认每个子句是否满足要求，再通过酶切等操作分解掉不满足当前时间步中的要求的片段，逐步筛选出满足所有约束条件的解。
  
DNA计算在解决SAT问题方面展现出独特优势：（1）大规模并行性。DNA计算的最大优势在于其天然的并行性。在一个试管环境中就可以同时处理万亿级别的DNA分子反应，相当于并行探索指数级数量的待选解。这种并行性使得原本需要指数时间的计算理论上可以在常数时间内完成。（2）超高存储密度。DNA分子具有较高的信息存储密度，1克DNA可以存储约215PB的数据。这对于存储大规模候选解非常有帮助。（3）生化操作的精确性。现代生物技术已经能够实现特异性DNA操作，包括精准的杂交、切割、复制，保证了计算过程的准确性。
  
虽然DNA计算也存在一些应用上的挑战，例如电泳速度较慢、提前准备DNA链较繁杂、需要溶液生化环境等，但DNA计算的信息编码、杂交匹配、并行性思想等核心思想启发了许多后续的研究工作和算法优化。
  
***
  
## 四 DNA计算解决SAT问题的经典实例
  
### 4.1 Qinghua Liu等(2000)：表面DNA计算解决4变量3-SAT问题
  
#### 实验简述
该实验首次在固体表面（金膜）上实现SAT计算，成功求解了一个包含4个变量、4个子句的3-SAT实例:<p align="center"><img src="https://latex.codecogs.com/gif.latex?\phi%20=%20(\bar{x_1}%20\vee%20\bar{x_2}%20\vee%20x_3)%20\wedge%20(x_1%20\vee%20\bar{x_2}%20\vee%20\bar{x_3})%20\wedge%20(x_1%20\vee%20x_2%20\vee%20\bar{x_3})%20\wedge%20(\bar{x_1}%20\vee%20x_2%20\vee%20x_4)"/></p>  
  
  
#### 实验原理   
1. 信息编码。合成2^4=16条不同的单链DNA分子，每条链的序列对应编码一个可能的变量赋值组合；合成4条互补的探针链，每条探针的序列存储一个子句的信息。如果某条DNA链满足某个子句，则该DNA链的序列与子句的序列互补。
  
<div align="center">
  <a href="https://imgchr.com/i/pZmDzHs">
    <img src="https://s41.ax1x.com/2025/12/06/pZmrZDJ.png" alt="pZmDzHs.png" width="50%" height="50%">/>
    <div style="color: #000; font-size: 0.9em; font-style: italic; margin-top: 1px;">
    图4.1.1：Liu等人(2000)实验的核心循环过程示意图
  </div>
  </a>
</div>
  
1. 循环筛选。如图4.1.1所示：针对问题中的每个子句，重复进行一套循环，模拟对所有候选解的以单个子句为步长的循环检查。（1）标记：加入探针，该探针与满足当前检测的子句的DNA序列互补，从而与表面上编码有“满足该子句赋值”信息的DNA链杂交结合，形成双链，而未能杂交的单链则代表不满足该子句。（2）销毁：加入核酸外切酶，用于降解单链DNA，在每一轮中销毁未被标记保护（不满足当前子句）的DNA链。（3）解标记：通过加热使双链DNA解链，移除探针，让表面剩余的DNA恢复为单链状态，以进行下一个子句的筛选。
2. 读取结果。完成所有子句循环后，只有包含满足所有子句正确答案的编码的DNA链会留在表面。使用荧光探针检测剩余序列，并将检测结果编码的答案信息还原，即可得到答案。
  
  
#### 实验流程 
实验流程共包含6个环节。如图4.1.2所示：（1）制备：将每种可能的赋值用一条15-mer单链核苷酸表示（即16种可能的信息链），并制备对应4个子句信息的4种探针链；（2）固定：将所有16条信息链固定在金表面形成DNA阵列；（3）杂交：依次加入可以与满足当前子句的DNA序列互补的探针，这些探针能“保护”符合当前子句赋值的单链；（4）分解：加入核酸外切酶Fok I，用于专一性降解单链DNA，在每一轮中销毁未被标记保护、即不满足当前轮次子句的DNA单链；（5）解标记：通过高温使得在（3）中吸附到信息链上的探针链脱离；（6）读取：经过数轮（3）到（5）的循环，只有包含满足所有子句正确答案的编码的DNA链会留在表面。使用荧光标记的读取探针检测剩余序列，通过PCA扩增剩余的DNA，然后即可确定其序列，还原出答案。
  
<div align="center">
  <a href="https://imgchr.com/i/pZmDzHs">
    <img src="https://s41.ax1x.com/2025/12/06/pZmDzHs.png" alt="pZmDzHs.png" width="55%" height="55%">/>
    <div style="color: #000; font-size: 0.9em; font-style: italic; margin-top: 1px;">
    图4.1.2：Liu等人(2000)的表面DNA计算实验流程示意图
  </div>
  </a>
</div>
  
  
  
关键技术细节：（1）使用硫醇修饰的DNA序列，使其能够牢固地附着在金膜表面；（2）每个变量的T/F状态通过特定的15-mer序列进行编码，例如“<img src="https://latex.codecogs.com/gif.latex?x_i"/>为真”用序列ATGCGATCGATCGAT表示；（3）探针设计采用了"保护性"策略，即只有满足特定子句条件的序列才能与探针形成稳定的双链结构，从而在后续酶切步骤中得到保护。
  
#### 实验结果  
成功在表面上唯一保留了满足所有子句的赋值，荧光信号清晰可辨，经过PCA扩增，正确答案对应DNA的浓度显著较高，成功鉴定出四个满足所有条件的正确答案，如图4.1.3。实验数据显示，正确与错误答案信号的信噪比在10到777之间，能够清晰地区分。
  
<div align="center">
  <a href="https://imgchr.com/i/pZmrnER">
    <img src="https://s41.ax1x.com/2025/12/06/pZmrnER.jpg" alt="pZmrnER.jpg"width="40%" height="40%" />
  </a>
  <div style="color: #000; font-size: 0.9em; font-style: italic; margin-top: 8px;">
    图4.1.3：表面DNA计算实验结果中的DNA浓度分布
  </div>
</div>
  
#### 实验意义与质疑
- 开创性意义：首次证明SAT可在固体表面以全自动、并行方式求解；提出了新范式，首次演示在固体表面通过可编程的生化操作解决计算问题，开创了“表面DNA计算”分支；展示了DNA计算的高度并行性，一次性操作海量分子（验证所有候选解），体现出强大的并行计算能力。
  
- 技术局限：实验中每个变量需要15碱基的DNA序列来编码，随着问题规模扩大，合成难度和耗费时间将急剧上升；后续研究指出了操作可能存在的误差性，“标记”步骤可能不完全，“销毁”步骤也可能不彻底，导致错误残留。因此有质疑者认为，随着问题规模（变量数）增大，操作误差会累积，可能导致可满足问题被误判为不可满足。
  
- 方法扩展的挑战：虽然该方法成功解决了4变量SAT问题，但将此方法扩展到更大规模问题将面临显著挑战。当变量数超过10时，实验复杂度会呈指数级增长，操作误差的累积效应将严重影响结果可靠性。
  
### 4.2 Ravinderjit S. Braich等(2002)：贴纸模型解决20变量3-SAT
  
#### 实验简述   
该实验采用“贴纸模型”（Sticker Model）和凝胶电泳等创新方法，成功求解了包含20个变量的3-SAT问题。（问题实例及答案见文末）
  
#### 实验原理  
  
本实验的核心原理基于“贴纸模型”，旨在通过可逆的核酸杂交过程来动态编码和操作信息。构建长链单链DNA“记忆链”，在其不同部位可与不同短链DNA“贴纸”结合，结合与否编码了对应的布尔信息。贴纸链通过碱基互补配对与记忆链结合，该过程可通过改变温度或化学环境实现可逆的“粘贴”与“擦除”。已经初始化的、包含所有可能解的DNA分子池，被引入一个由多个“分离器”模块构成的网络中，如图4.2.1。每个分离器模块包含某一特定贴纸DNA，设计为仅捕获满足某一特定子句的DNA分子以模拟子句逻辑门，通过连续的“过滤”操作完成对所有约束条件的并行性验证。最终留下的满足所有约束条件的DNA分子，即包含了所求解的信息。
  
<div style="text-align: center;">
  <a href="https://imgchr.com/i/pZmyU9P">
    <img src="https://s41.ax1x.com/2025/12/06/pZmyU9P.png" alt="pZmyU9P.png" style="width: 65%; height: 65%; margin: auto;" />
  </a>
  <div style="color: #000; font-size: 0.9em; font-style: italic; margin-top: 8px;">
    图4.2.1：Braich等人(2002)的贴纸模型中的分离器
  </div>
</div>
  
#### 实验流程  
1. 库的初始化：合成经过设计的长度固定的（300 nt）记忆链库。随后根据所求解的SAT问题具体形式，通过一系列杂交步骤，将对应的贴纸链粘贴到记忆链上，从而构建出一个代表所有2²⁰种可能赋值组合的初始候选解空间。
2. 分离器模块的设计：制备一系列填充有聚丙烯酰胺凝胶的玻璃模块，每个模块的凝胶中固定有特定DNA探针序列，该序列编码了待求解的3-SAT问题中某一个子句的满足条件，会与能与满足此子句所对应的记忆链区域互补。
3. 循环筛选计算：如图4.2.1，初始库被加载到第一个分离器模块中，随后开始如下循环过程：（1）电泳驱动与特异性捕获：在电场作用下，DNA分子在模块间定向流动。当流经模块时，序列能与凝胶内固定探针互补的DNA链被捕获并滞留在凝胶中；不满足的链则顺利通过，进入废液池。（2）洗脱与转移：完成一个模块的筛选后，通过加热使被捕获的DNA链从探针上脱离。随后，剩下的经过筛选的液体，即满足上一轮子句的解的集合，被电泳驱动至下一个模块，接受下一轮筛选。
4. 输出检测与分析：经过所有模块筛选后，最终流出的分子仅包含满足全部子句的正确解。通过PRC反应对微量的最终产物进行扩增，随后解码其对应的变量，得到满足所有子句的解。
  
#### 实验结果 
经过筛选后，最终的分子产物编码的变量赋值完全满足给定的20变量3-SAT问题的所有约束条件，如图4.2.2所示。实验成功地从1048576种可能性中找到了唯一解。虽然过程中存在因非完全捕获或洗脱导致的分子脱落，但通过精心的编码设计与流程优化，有效信号与噪声背景的信噪比足以证明实验的有效性。
  
<div align="center">
  <a href="https://imgchr.com/i/pZm69UA">
    <img src="https://s41.ax1x.com/2025/12/06/pZm69UA.png" alt="pZm69UA.png" style="width: 80%;" />
  </a>
  <div style="color: #000; font-size: 0.9em; font-style: italic; margin-top: 8px;">
    图4.2.2：Braich等人(2002)的贴纸模型中最终产物对应编码的信号
  </div>
</div>
  
#### 实验意义与质疑  
- 意义
  - 该实验是首次使用非硅基方式解决变量数达到20的NP完全问题，极大地提升了DNA计算处理问题的规模，证明了用分子系统执行复杂多步离散算法的可行性。
  - 实验提出的“贴纸模型”等架构为DNA计算提供了一种模块化、可编程、自动化的范式，扩展了DNA计算的通用性。
  - 该工作极大鼓舞并引领了后续十余年的DNA计算研究，激发了大量解决其他NP问题（如图着色、集合覆盖等）的算法设计。
- 质疑
  - 尽管解决了20变量问题，进一步扩大规模仍面临指数级增长的挑战。每一步的固有生化误差会逐级累积，最终导致信噪比过低，结果不可靠。后续分析指出，此类方法在应对更大规模NP完全问题时存在原理性限制。
  - 整个计算过程耗时数小时至数天，其计算时间主要消耗在缓慢的电泳迁移和温控循环上，时间效率远低于电子计算机对于同类问题的求解。
  - 针对特定问题设计、合成和验证的记忆链和探针序列，本身是一项繁琐的任务。复杂的编码方案也增加了设计的难度和复杂性。
  
***
  
## 五 DNA计算的开源仿真代码及添加可视化
  
  
<table style="width: 100%; border: none;">
<tr>
<td style="width: 0%; vertical-align: top; padding-right: 0px;">
  
<div align="center">
  <a href="https://imgchr.com/i/pZ13LCT">
    <img src="https://s41.ax1x.com/2025/12/18/pZ13LCT.png" alt="pZm69UA.png" style="width: 80%;" />
  </a>
  <div style="color: #000; font-size: 1.0em; font-style: italic; margin-top: 8px;">
    图5.1:原开源代码模拟DNA计算解决SAT问题流程
  </div>
</div>
  
  
</td>
<td style="width: 10%; vertical-align: top;">
  
### 5.1 Python仿真框架[^3]
  
该Python仿真框架用于模拟DNA计算中的关键生化操作。采用面向对象设计，包含DNA链的初始化、扩增、分离、筛选等操作。仿真流程如下：
  
1. 初始化：创建初始DNA链集合，代表所有可能的变量赋值组合（候选解）
2. 生化操作：依次执行PCR扩增、分离、筛选等操作
3. 结果验证：通过凝胶电泳和序列分析验证最终结果是否正确
  
主要功能：
- 碱基操作：互补序列的计算
- 扩增：模拟PCR扩增，增多特定序列
- 分离：按长度等特征分离DNA链
- 筛选：基于特定条件筛选DNA链
- 连接：模拟DNA连接酶作用
- 电泳：模拟凝胶电泳分离过程
  
**Main中代码执行流程**：
  
```python
# 创建DNA计算模拟器实例
calc = DNACalculator()
# 定义初始DNA链集合，代表不同的计算输入
initial_strands = pass
# 初始化DNA链集合
calc.initialize_strands(initial_strands)
# 扩增特定序列，模拟增加特定输入的权重
calc.amplify("ABCDEFGHIJKLMNOPQRSTUVWXYZ")
# 按长度分离DNA链
length_groups = calc.separate_by_length()
# 选择特定长度的DNA链
selected_by_length = calc.select_by_length(26)
# 根据模式分离DNA链
with_pattern, without_pattern = calc.separate_by_pattern("ABC")
# 选择包含特定模式的DNA链
selected_by_pattern = calc.select_by_pattern("ABC")
# 连接两条DNA链，模拟逻辑运算
ligated = calc.ligate("ABC", "DEF")
# 进行PCR扩增，模拟信号放大
calc.pcr("ABC", "DEF", cycles=3)
# 进行凝胶电泳，选择特定长度的DNA链
calc.gel_electrophoresis(26)
# 打印整个计算过程的所有中间结果
calc.print_intermediate_results()
# 可视化DNA链演化过程
calc.visualize_strand_evolution()
```
</td>
</tr>
</table>
  
### 5.2 基本类设计
  
核心类为`DNACalculator`，封装了DNA链的初始化、扩增、分离、选择等操作。
```python
    class DNACalculator:
        def __init__(self):
            # 存储当前所有的DNA序列
            self.strands = []
            # 存储中间结果
            self.intermediate_results = []
```
  
### 5.3 底层生化操作
  
#### 5.3.1 互补序列与反向互补序列
  
仿真器实现了基于碱基互补配对原则的方法`complement()`和`reverse_complement()`：
  
```python
    # 定义DNA碱基映射关系
    base_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    def complement(sequence):
        """获取互补序列"""
        return ''.join([base_map.get(base, base) for base in sequence])
    def reverse_complement(sequence):
        """获取反向互补序列"""
        return complement(sequence)[::-1]
```
  
#### 5.3.2 PCR扩增仿真
  
PCR是DNA计算中的关键扩增技术，通过`amplify()`方法实现：
  
```python
    def amplify(self, target_sequence, fold=2):
        """
        PCR扩增特定序列
        参数:
            target_sequence (str): 需要扩增的目标序列
            fold (int): 扩增倍数，默认为2
        """
        count = self.strands.count(target_sequence)
        self.strands.extend([target_sequence] * count * fold)
        self.intermediate_results.append((f'amplify {target_sequence} {fold}x', self.strands.copy()))
```
  
#### 5.3.3 分离与筛选操作
  
提供根据长度或模式分离DNA链的方法：
  
```python
    def separate_by_length(self):
        """按长度分离DNA链"""
        length_groups = {}
        for strand in self.strands:
            length = len(strand)
            if length not in length_groups:
                length_groups[length] = []
            length_groups[length].append(strand)
        return length_groups
```
```python
    def separate_by_pattern(self, pattern):
        """根据模式分离DNA链"""
        with_pattern = []
        without_pattern = []
        for strand in self.strands:
            if pattern in strand:
                with_pattern.append(strand)
            else:
                without_pattern.append(strand)
        return with_pattern, without_pattern
```
  
### 5.4 高级生化操作
  
#### 5.4.1 DNA连接
  
DNA连接用于组合不同的DNA片段，以创建新的基因序列：
  
```python
    def ligate(self, sequence1, sequence2):
        """连接两条DNA片段"""
        ligated = sequence1 + sequence2
        self.strands.append(ligated)
        self.intermediate_results.append((f'ligate {sequence1} + {sequence2}', ligated))
        return ligated
```
  
#### 5.4.2 凝胶电泳
  
通过凝胶电泳，按长度分离DNA分子：
  
```python
    def gel_electrophoresis(self, length):
        """模拟凝胶电泳，筛选特定长度的DNA片段"""
        self.strands = [strand for strand in self.strands if len(strand) == length]
        self.intermediate_results.append((f'gel electrophoresis select length {length}', self.strands.copy()))
        return self.strands
```
  
#### 5.4.3 PCR完整过程
  
完整的PCR过程，包括引物结合和片段扩增：
  
```python
    def pcr(self, primer1, primer2, cycles=10):
        """模拟PCR完整过程"""
        reverse_primer1 = reverse_complement(primer1)
        reverse_primer2 = reverse_complement(primer2)
  
        for cycle in range(cycles):
            new_strands = []
            for strand in self.strands:
                if primer1 in strand and reverse_primer2 in strand:
                    start_index = strand.find(primer1)
                    end_index = strand.find(reverse_primer2) + len(reverse_primer2)
                    if start_index < end_index:
                        amplified_segment = strand[start_index:end_index]
                        new_strands.append(amplified_segment)
            self.strands.extend(new_strands)
```
  
### 5.5 过程结果追踪
  
完整的追踪方法，便于分析过程：
  
```python
    def print_strands(self):
        """打印当前所有DNA链及信息"""
        print("Current strands:")
        strand_counts = Counter(self.strands)
        for strand, count in strand_counts.items():
            print(f"  {strand}: {count}")
        print()
    ```
    ```python
    def print_intermediate_results(self):
        """打印所有中间结果"""
        print("Intermediate results:")
        for step, result in self.intermediate_results:
            print(f"  Step: {step}")
            pass
```
  
### 5.6 添加可视化
  
为原代码添加两个可视化：DNA总数及不同种类DNA链的数量变化，用于清晰地展示不同步骤中DNA链的演化过程。
  
```python
  
def visualize_strand_evolution(self):
        """
        可视化DNA链的演化过程
        """
        if not self.intermediate_results:
            print("No intermediate results to visualize")
            return
  
        steps = []
        step_names = []
        strand_counts = []
        strand_distributions = []  # 存储每个步骤中不同DNA序列的数量
  
        for i, (step_name, result) in enumerate(self.intermediate_results):
            steps.append(i)
            step_names.append(step_name)
  
            if isinstance(result, list):
                strand_counts.append(len(result))
                # 统计每种DNA序列的总数
                counter = Counter(result)
                strand_distributions.append(counter)
            elif isinstance(result, dict):
                # 如果是字典格式，就计算所有链的总数
                total_count = 0
                all_strands = []
                for key, value in result.items():
                    if isinstance(value, list):
                        total_count += len(value)
                        all_strands.extend(value)
                strand_counts.append(total_count)
                # 统计每种DNA序列的总数
                counter = Counter(all_strands)
                strand_distributions.append(counter)
            else:
                strand_counts.append(0)
                strand_distributions.append(Counter())
  
```
  
```python
        # 创建图表
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 12))
  
        # 图表1: DNA链总数变化（折线图）
        x = np.arange(len(steps))
  
        ax1.plot(x, strand_counts, marker='o', linestyle='-', linewidth=2, 
                label='Total DNA Strands', color='skyblue', markersize=6)
  
        ax1.set_xlabel('Processing Steps')
        ax1.set_ylabel('Number of DNA Strands')
        ax1.set_title('Total DNA Strands Evolution')
        ax1.set_xticks(x)
        ax1.set_xticklabels(step_names, rotation=45, ha='right')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
  
        # 在折线图上显示具体数值
        for i, v in enumerate(strand_counts):
            ax1.text(i, v + 0.5, str(v), ha='center', va='bottom', fontsize=9)
  
```
```python
        # 图表2: 每种DNA链的个数（堆叠柱状图）
        # 获取所有步骤中出现过的所有DNA序列
        all_sequences = set()
        for counter in strand_distributions:
            all_sequences.update(counter.keys())
        all_sequences = sorted(list(all_sequences))
  
        # 为每种DNA序列分配颜色
        colors = plt.cm.Set3(np.linspace(0, 1, len(all_sequences)))
  
        # 准备堆叠数据
        bottom = np.zeros(len(steps))
        for seq_idx, sequence in enumerate(all_sequences):
            counts = []
            for counter in strand_distributions:
                counts.append(counter.get(sequence, 0))
  
            ax2.bar(step_names, counts, bottom=bottom, label=sequence, 
                   color=colors[seq_idx % len(colors)], alpha=0.8)
            bottom += np.array(counts)
  
        ax2.set_xlabel('Processing Steps')
        ax2.set_ylabel('Number of DNA Strands')
        ax2.set_title('DNA Strand Distribution by Type')
        ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        ax2.set_xticklabels(step_names, rotation=45, ha='right')
        ax2.grid(True, alpha=0.3)
  
        plt.tight_layout()
        plt.show()
        plt.close()
```
  
### 5.7 仿真示例与结果展示
  
一个完整的DNA计算仿真示例，展示了从初始化到最终结果的全过程及对应可视化：
  
```python
    def main():
        calc = DNACalculator()
  
        # 初始化DNA链集合
        initial_strands = [
            "ABCDEFGHIJKLMNOPQRSTUVWXYZ",      # 输入A
            "XYZABCDEFGHIJKLMNOPQRSTUVW",      # 输入B
            "ABCDEFGHIJKLMNOPQRSTUVWXY",       # 输入C
            "QRSTUVWXYZABCDEFGHIJKLMNOP"       # 输入D
        ]
  
        calc.initialize_strands(initial_strands)
        calc.print_strands()
```
```python
  
        calc.amplify("ABCDEFGHIJKLMNOPQRSTUVWXYZ")
        calc.separate_by_length()
        calc.select_by_length(26)
        calc.separate_by_pattern("ABC")
        calc.ligate("ABC", "DEF")
        calc.pcr("ABC", "DEF", cycles=3)
        calc.gel_electrophoresis(26)
  
        calc.print_intermediate_results()
        #完整输出结果详见文末
```
  
<div align="center">
  <a href="https://imgchr.com/i/pZ1RwRS">
    <img src="https://s41.ax1x.com/2025/12/19/pZ1RwRS.png" alt="pZm69UA.png" style="width: 100%;" />
  </a>
  <div style="color: #000; font-size: 1.0em; font-style: italic; margin-top: 8px;">
    图5.7:原开源代码模拟DNA计算解决SAT问题流程
  </div>
</div>
  
### 5.8 DNA计算仿真的意义
  
Python仿真DNA计算具有如下重要意义。（1）支持复杂SAT问题求解，可以利用计算机处理多变量SAT问题；（2）可视化计算过程：利用matplotlib可视化清晰展示DNA链的演化过程与模拟溶液环境；（3）性能优化：针对大规模问题优化算法效率较高，不必制备繁杂的DNA分子与引物；（4）环境稳定：不需考虑溶液生化环境，直接在模拟中计算，方便快捷。
  
***
  
## 六 DNA计算的后续发展与启发意义
  
自2002年后，DNA计算在SAT及更广义的约束满足问题领域持续发展，主要方向包括：
  
1. 自动化与微流控：DNA逻辑门电路 → 完全自动化的20变量SAT求解器
  
2. 分子编程语言与编译器：CRN-to-DNA编译器
  
3. 与其他新兴计算融合：如DNA+遗传算法、图神经网络（GNN）
  
4. 规模化尝试：利用DNA origami解决了一个70变量的随机3-SAT实例 *（施一公，2021）*
  
尽管目前DNA计算在实际求解大规模SAT实例时仍无法与顶级求解器竞争，但其理论意义深远：证明了生物分子可实现复杂逻辑计算；大规模并行性为破解指数时间限制提供了新思路；启发了量子计算中“全局搜索”机制的设计；为生物学与计算学科交叉提供了范式；启发了作者对利用元胞自动机、遗传算法解决逻辑运算的思考与对离散数学学科的兴趣。（是真的）
  
***
  
## 七 参考文献及相关材料
Liu, Q., Wang, L., Frutos, A. G., Condon, A. E., Corn, R. M., & Smith, L. M. (2000). DNA computing on surfaces. *Nature*, 403(6766), 175–179. https://doi.org/10.1038/35003155
  
Braich, R. S., Chelyapov, N., Johnson, C., Rothemund, P. W., & Adleman, L. M. (2002). Solution of a 20-variable 3-SAT problem on a DNA computer. *Science*, 296(5567), 499–502. https://doi.org/10.1126/science.1069528
  
Adleman, L. M. (1994). Molecular computation of solutions to combinatorial problems. Science, *266*(5187), 1021–1024. https://doi.org/10.1126/science.7973651
  
Zobront, J. (n.d.). DNA computing simulator [Computer software]. GitHub. Retrieved from https://github.com/zobront/dna-computing-simulator
  
  
[^1]:https://doi.org/10.1038/35003155
[^2]:https://doi.org/10.1126/science.1069528
[^3]:https://github.com/zobront/dna-computing-simulator
  
<div style="text-align: left;">
<div style="color: #000; font-size: 1.0em; font-style: italic; margin-top: 8px;">
    图1:20SAT问题中的问题与解
  </div>
  <a href="https://imgchr.com/i/pZms0oR">
    <img src="https://s41.ax1x.com/2025/12/06/pZms0oR.md.png" alt="pZms0oR.md.png" style="width: 60%;" />
  </a>
  
</div>
  
代码块1：Python代码仿真后输出结果
  
```
Initializing strands...
Current strands:
  ABCDEFGHIJKLMNOPQRSTUVWXYZ: 1
  XYZABCDEFGHIJKLMNOPQRSTUVW: 1
  ABCDEFGHIJKLMNOPQRSTUVWXY: 1
  QRSTUVWXYZABCDEFGHIJKLMNOP: 1
  ```
  ```
Amplifying specific strands...
Current strands:
  ABCDEFGHIJKLMNOPQRSTUVWXYZ: 3
  XYZABCDEFGHIJKLMNOPQRSTUVW: 1
  ABCDEFGHIJKLMNOPQRSTUVWXY: 1
  QRSTUVWXYZABCDEFGHIJKLMNOP: 1
  ```
  ```
Separating strands by length...
{26: ['ABCDEFGHIJKLMNOPQRSTUVWXYZ', 
'XYZABCDEFGHIJKLMNOPQRSTUVW',
'QRSTUVWXYZABCDEFGHIJKLMNOP', 
'ABCDEFGHIJKLMNOPQRSTUVWXYZ', 
'ABCDEFGHIJKLMNOPQRSTUVWXYZ'], 25: 
['ABCDEFGHIJKLMNOPQRSTUVWXY']}
```
```
Selecting strands of specific length...
['ABCDEFGHIJKLMNOPQRSTUVWXYZ',
'XYZABCDEFGHIJKLMNOPQRSTUVW',
'QRSTUVWXYZABCDEFGHIJKLMNOP',
'ABCDEFGHIJKLMNOPQRSTUVWXYZ',
'ABCDEFGHIJKLMNOPQRSTUVWXYZ']
```
```
Separating strands by pattern...
With pattern 'ABC': 
['ABCDEFGHIJKLMNOPQRSTUVWXYZ',
'XYZABCDEFGHIJKLMNOPQRSTUVW',
'ABCDEFGHIJKLMNOPQRSTUVWXY',
'QRSTUVWXYZABCDEFGHIJKLMNOP',
'ABCDEFGHIJKLMNOPQRSTUVWXYZ',
'ABCDEFGHIJKLMNOPQRSTUVWXYZ']
Without pattern 'ABC': []
```
```
Selecting strands with specific pattern...
['ABCDEFGHIJKLMNOPQRSTUVWXYZ',
'XYZABCDEFGHIJKLMNOPQRSTUVW',
'ABCDEFGHIJKLMNOPQRSTUVWXY',
'QRSTUVWXYZABCDEFGHIJKLMNOP',
'ABCDEFGHIJKLMNOPQRSTUVWXYZ',
'ABCDEFGHIJKLMNOPQRSTUVWXYZ']
```
```
Ligating two strands...
Current strands:
  ABCDEFGHIJKLMNOPQRSTUVWXYZ: 3
  XYZABCDEFGHIJKLMNOPQRSTUVW: 1
  ABCDEFGHIJKLMNOPQRSTUVWXY: 1
  QRSTUVWXYZABCDEFGHIJKLMNOP: 1
  ABCDEF: 1
  ```
  ```
Performing PCR...
Current strands:
  ABCDEFGHIJKLMNOPQRSTUVWXYZ: 3
  XYZABCDEFGHIJKLMNOPQRSTUVW: 1
  ABCDEFGHIJKLMNOPQRSTUVWXY: 1
  QRSTUVWXYZABCDEFGHIJKLMNOP: 1
  ABCDEF: 1
  ```
```
Performing gel electrophoresis...
Current strands:
  ABCDEFGHIJKLMNOPQRSTUVWXYZ: 3
  XYZABCDEFGHIJKLMNOPQRSTUVW: 1
  QRSTUVWXYZABCDEFGHIJKLMNOP: 1
  ```
```
Intermediate results:
  Step: initialize
    ['ABCDEFGHIJKLMNOPQRSTUVWXYZ',
    'XYZABCDEFGHIJKLMNOPQRSTUVW',
    'ABCDEFGHIJKLMNOPQRSTUVWXY',
    'QRSTUVWXYZABCDEFGHIJKLMNOP']
  ```
  ```
  Step: amplify ABCDEFGHIJKLMNOPQRSTUVWXYZ 2x
    ['ABCDEFGHIJKLMNOPQRSTUVWXYZ',
    'XYZABCDEFGHIJKLMNOPQRSTUVW',
    'ABCDEFGHIJKLMNOPQRSTUVWXY',
    'QRSTUVWXYZABCDEFGHIJKLMNOP',
    'ABCDEFGHIJKLMNOPQRSTUVWXYZ',
    'ABCDEFGHIJKLMNOPQRSTUVWXYZ']
  ```
  ```
  Step: separate by length
    26: ['ABCDEFGHIJKLMNOPQRSTUVWXYZ',
    'XYZABCDEFGHIJKLMNOPQRSTUVW',
    'QRSTUVWXYZABCDEFGHIJKLMNOP',
    'ABCDEFGHIJKLMNOPQRSTUVWXYZ',
    'ABCDEFGHIJKLMNOPQRSTUVWXYZ']
    25: ['ABCDEFGHIJKLMNOPQRSTUVWXY']
  ```
  ```
  Step: select length 26
    ['ABCDEFGHIJKLMNOPQRSTUVWXYZ',
    'XYZABCDEFGHIJKLMNOPQRSTUVW',
    'QRSTUVWXYZABCDEFGHIJKLMNOP',
    'ABCDEFGHIJKLMNOPQRSTUVWXYZ',
    'ABCDEFGHIJKLMNOPQRSTUVWXYZ']
  ```
  ```
  Step: separate by pattern ABC
    with: ['ABCDEFGHIJKLMNOPQRSTUVWXYZ',
    'XYZABCDEFGHIJKLMNOPQRSTUVW',
    'ABCDEFGHIJKLMNOPQRSTUVWXY',
    'QRSTUVWXYZABCDEFGHIJKLMNOP',
    'ABCDEFGHIJKLMNOPQRSTUVWXYZ',
    'ABCDEFGHIJKLMNOPQRSTUVWXYZ']
    without: []
  ```
  ```
  Step: select pattern ABC
    ['ABCDEFGHIJKLMNOPQRSTUVWXYZ',
    'XYZABCDEFGHIJKLMNOPQRSTUVW',
    'ABCDEFGHIJKLMNOPQRSTUVWXY',
    'QRSTUVWXYZABCDEFGHIJKLMNOP',
    'ABCDEFGHIJKLMNOPQRSTUVWXYZ',
    'ABCDEFGHIJKLMNOPQRSTUVWXYZ']
  ```
  ```
  Step: ligate ABC + DEF
  Step: pcr with primers ABC, DEF for 3 cycles
    ['ABCDEFGHIJKLMNOPQRSTUVWXYZ',
    'XYZABCDEFGHIJKLMNOPQRSTUVW',
    'ABCDEFGHIJKLMNOPQRSTUVWXY',
    'QRSTUVWXYZABCDEFGHIJKLMNOP',
    'ABCDEFGHIJKLMNOPQRSTUVWXYZ',
    'ABCDEFGHIJKLMNOPQRSTUVWXYZ', 'ABCDEF']
  ```
  ```
  Step: gel electrophoresis select length 26
    ['ABCDEFGHIJKLMNOPQRSTUVWXYZ',
    'XYZABCDEFGHIJKLMNOPQRSTUVW',
    'QRSTUVWXYZABCDEFGHIJKLMNOP',
    'ABCDEFGHIJKLMNOPQRSTUVWXYZ',
    'ABCDEFGHIJKLMNOPQRSTUVWXYZ']
  ```
  