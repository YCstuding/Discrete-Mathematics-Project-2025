import random
from collections import Counter
import matplotlib.pyplot as plt
import numpy as np

# 定义DNA碱基映射关系，用于获取互补碱基
# 标准DNA碱基包括 A(腺嘌呤), T(胸腺嘧啶), C(胞嘧啶), G(鸟嘌呤)
base_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

def complement(sequence):
    """
    获取DNA序列的互补序列
    
    对于非标准DNA碱基，保持原样不变
    
    参数:
        sequence (str): 原始DNA序列
        
    返回:
        str: 互补DNA序列
    """
    # 对于每个碱基，如果在映射表中则替换为互补碱基，否则保持原样
    return ''.join([base_map.get(base, base) for base in sequence])

def reverse_complement(sequence):
    """
    获取DNA序列的反向互补序列
    
    参数:
        sequence (str): 原始DNA序列
        
    返回:
        str: 反向互补DNA序列
    """
    return complement(sequence)[::-1]

class DNACalculator:
    """
    DNA计算模拟器类
    模拟基于DNA链的计算过程，包括各种生物化学操作
    """
    
    def __init__(self):
        """
        初始化DNA计算模拟器
        """
        # 存储当前所有的DNA序列
        self.strands = []
        # 存储中间结果，用于调试和追踪
        self.intermediate_results = []

    def initialize_strands(self, initial_strands):
        """
        初始化DNA链集合
        
        参数:
            initial_strands (list): 初始DNA序列列表
        """
        self.strands = initial_strands.copy()
        # 记录初始状态
        self.intermediate_results.append(('initialize', self.strands.copy()))

    def amplify(self, target_sequence, fold=2):
        """
        扩增特定DNA序列（模拟PCR扩增）
        
        参数:
            target_sequence (str): 需要扩增的目标序列
            fold (int): 扩增倍数，默认为2
        """
        # 查找目标序列在当前链中的数量
        count = self.strands.count(target_sequence)
        # 将目标序列添加fold*count次到链集合中
        self.strands.extend([target_sequence] * count * fold)
        # 记录扩增操作结果
        self.intermediate_results.append((f'amplify {target_sequence} {fold}x', self.strands.copy()))

    def separate_by_length(self):
        """
        按长度分离DNA链
        将不同长度的DNA链分组存储
        
        返回:
            dict: 以长度为键，DNA链列表为值的字典
        """
        # 使用collections.Counter统计各长度序列的数量
        length_groups = {}
        for strand in self.strands:
            length = len(strand)
            if length not in length_groups:
                length_groups[length] = []
            length_groups[length].append(strand)
        
        # 记录分离结果
        self.intermediate_results.append(('separate by length', length_groups))
        return length_groups

    def separate_by_pattern(self, pattern):
        """
        根据模式分离DNA链
        分离包含和不包含特定模式的DNA链
        
        参数:
            pattern (str): 分离依据的模式序列
            
        返回:
            tuple: (包含模式的链列表, 不包含模式的链列表)
        """
        with_pattern = []
        without_pattern = []
        # 遍历所有链，根据是否包含pattern进行分类
        for strand in self.strands:
            if pattern in strand:
                with_pattern.append(strand)
            else:
                without_pattern.append(strand)
        
        # 记录分离结果
        self.intermediate_results.append((f'separate by pattern {pattern}', 
                                        {'with': with_pattern, 'without': without_pattern}))
        return with_pattern, without_pattern

    def select_by_length(self, length):
        """
        选择特定长度的DNA链
        
        参数:
            length (int): 目标长度
            
        返回:
            list: 指定长度的所有DNA链
        """
        # 筛选出指定长度的链
        selected = [strand for strand in self.strands if len(strand) == length]
        # 记录选择结果
        self.intermediate_results.append((f'select length {length}', selected))
        return selected

    def select_by_pattern(self, pattern):
        """
        选择包含特定模式的DNA链
        
        参数:
            pattern (str): 目标模式
            
        返回:
            list: 包含指定模式的所有DNA链
        """
        # 筛选出包含指定模式的链
        selected = [strand for strand in self.strands if pattern in strand]
        # 记录选择结果
        self.intermediate_results.append((f'select pattern {pattern}', selected))
        return selected

    def ligate(self, sequence1, sequence2):
        """
        连接两条DNA链（模拟DNA连接酶作用）
        
        参数:
            sequence1 (str): 第一条DNA链
            sequence2 (str): 第二条DNA链
            
        返回:
            str: 连接后的DNA链
        """
        ligated = sequence1 + sequence2
        # 将连接产物加入链集合
        self.strands.append(ligated)
        # 记录连接操作结果
        self.intermediate_results.append((f'ligate {sequence1} + {sequence2}', ligated))
        return ligated

    def gel_electrophoresis(self, length):
        """
        模拟凝胶电泳，筛选特定长度的DNA链
        
        参数:
            length (int): 目标长度
            
        返回:
            list: 电泳后保留的DNA链
        """
        # 保留指定长度的链，移除其他长度的链
        self.strands = [strand for strand in self.strands if len(strand) == length]
        # 记录电泳结果
        self.intermediate_results.append((f'gel electrophoresis select length {length}', self.strands.copy()))
        return self.strands

    def pcr(self, primer1, primer2, cycles=10):
        """
        模拟聚合酶链反应(PCR)过程
        
        参数:
            primer1 (str): 正向引物
            primer2 (str): 反向引物
            cycles (int): PCR循环次数
        """
        # 获取引物的反向互补序列
        reverse_primer1 = reverse_complement(primer1)
        reverse_primer2 = reverse_complement(primer2)
        
        # 进行指定次数的PCR循环
        for cycle in range(cycles):
            new_strands = []
            # 遍历当前所有链
            for strand in self.strands:
                # 如果链同时包含两个引物的结合位点
                if primer1 in strand and reverse_primer2 in strand:
                    # 提取两个引物之间的片段并进行扩增
                    start_index = strand.find(primer1)
                    end_index = strand.find(reverse_primer2) + len(reverse_primer2)
                    if start_index < end_index:
                        amplified_segment = strand[start_index:end_index]
                        new_strands.append(amplified_segment)
            
            # 将新扩增的链加入集合
            self.strands.extend(new_strands)
        
        # 记录PCR结果
        self.intermediate_results.append((f'pcr with primers {primer1}, {primer2} for {cycles} cycles', 
                                         self.strands.copy()))

    def print_strands(self):
        """
        打印当前所有DNA链及其统计信息
        """
        print("Current strands:")
        # 统计各序列的数量
        strand_counts = Counter(self.strands)
        # 按序列显示及其数量
        for strand, count in strand_counts.items():
            print(f"  {strand}: {count}")
        print()

    def print_intermediate_results(self):
        """
        打印所有中间结果，用于追踪计算过程
        """
        print("Intermediate results:")
        for step, result in self.intermediate_results:
            print(f"  Step: {step}")
            if isinstance(result, dict):
                # 如果结果是字典格式，分别打印各项
                for key, value in result.items():
                    if isinstance(value, list):
                        print(f"    {key}: {value}")
                    else:
                        print(f"    {key}: {value}")
            elif isinstance(result, list):
                # 如果结果是列表格式，直接打印
                print(f"    {result}")
            else:
                # 其他格式直接打印
                print(f"    {result}")
            print()

    def visualize_strand_evolution(self):
        """
        可视化DNA链的演化过程
        """
        if not self.intermediate_results:
            print("No intermediate results to visualize")
            return
        
        # 准备数据
        steps = []
        step_names = []
        strand_counts = []
        strand_distributions = []  # 存储每个步骤中不同DNA序列的分布
        
        for i, (step_name, result) in enumerate(self.intermediate_results):
            steps.append(i)
            step_names.append(step_name)
            
            if isinstance(result, list):
                strand_counts.append(len(result))
                # 统计每种DNA序列的个数
                counter = Counter(result)
                strand_distributions.append(counter)
            elif isinstance(result, dict):
                # 如果是字典格式，计算所有链的总数
                total_count = 0
                all_strands = []
                for key, value in result.items():
                    if isinstance(value, list):
                        total_count += len(value)
                        all_strands.extend(value)
                strand_counts.append(total_count)
                # 统计每种DNA序列的个数
                counter = Counter(all_strands)
                strand_distributions.append(counter)
            else:
                strand_counts.append(0)
                strand_distributions.append(Counter())
        
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




def main():
    """
    主函数，演示DNA计算模拟器的使用流程
    模拟解决一个简单的计算问题
    """
    # 创建DNA计算模拟器实例
    calc = DNACalculator()
    
    # 定义初始DNA链集合，代表不同的计算输入
    # 每个字母序列代表一个数值或状态
    initial_strands = [
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ",      # 输入A
        "XYZABCDEFGHIJKLMNOPQRSTUVW",      # 输入B
        "ABCDEFGHIJKLMNOPQRSTUVWXY",        # 输入C
        "QRSTUVWXYZABCDEFGHIJKLMNOP"       # 输入D
    ]
    
    print("Initializing strands...")
    # 初始化DNA链集合
    calc.initialize_strands(initial_strands)
    calc.print_strands()
    
    print("Amplifying specific strands...")
    # 扩增特定序列，模拟增加特定输入的权重
    calc.amplify("ABCDEFGHIJKLMNOPQRSTUVWXYZ")
    calc.print_strands()
    
    print("Separating strands by length...")
    # 按长度分离DNA链
    length_groups = calc.separate_by_length()
    print(length_groups)
    print()
    
    print("Selecting strands of specific length...")
    # 选择特定长度(26)的DNA链
    selected_by_length = calc.select_by_length(26)
    print(selected_by_length)
    print()
    
    print("Separating strands by pattern...")
    # 根据模式"ABC"分离DNA链
    with_pattern, without_pattern = calc.separate_by_pattern("ABC")
    print("With pattern 'ABC':", with_pattern)
    print("Without pattern 'ABC':", without_pattern)
    print()
    
    print("Selecting strands with specific pattern...")
    # 选择包含模式"ABC"的DNA链
    selected_by_pattern = calc.select_by_pattern("ABC")
    print(selected_by_pattern)
    print()
    
    print("Ligating two strands...")
    # 连接两条DNA链，模拟逻辑运算
    ligated = calc.ligate("ABC", "DEF")
    calc.print_strands()
    
    print("Performing PCR...")
    # 进行PCR扩增，模拟信号放大
    # 使用"ABC"作为正向引物，"DEF"的反向互补作为反向引物
    calc.pcr("ABC", "DEF", cycles=3)
    calc.print_strands()
    
    print("Performing gel electrophoresis...")
    # 进行凝胶电泳，选择特定长度(26)的DNA链
    calc.gel_electrophoresis(26)
    calc.print_strands()
    
    # 打印整个计算过程的所有中间结果
    calc.print_intermediate_results()
    
    # 可视化DNA链演化过程
    calc.visualize_strand_evolution()


# 程序入口点
if __name__ == "__main__":
    main()