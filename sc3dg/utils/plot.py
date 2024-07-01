import os
import matplotlib.pyplot as plt
import numpy as np

def DrawCellInternerCount(CellCounts = None, ChromLenList = None,save=''):
    # Fig.1 Cells Chromatain Interner Count and Proporation
    ###Get Chrom Length
    values_list = [value for value in ChromLenList[0].values()]
    dataPoints = sum(np.multiply(values_list, values_list))
    ###Get propotion
    CellCounts = sorted(CellCounts)
    CellProporty = CellCounts / dataPoints
    ###Log10 count
    CellCounts = np.log10(np.array(CellCounts) + 1)
    ###Draw Plot
    x = list(range(1, len(CellCounts)+1))
    # 创建图形和坐标轴对象
    fig, ax1 = plt.subplots()
    # 绘制第一个列表的数据
    ax1.plot(x, CellCounts, 'ro-', label='Count')
    ax1.set_xlabel('Cell Number')
    ax1.set_ylabel('1*10^n interactions')
    ax1.tick_params(axis='y')
    # 创建第二个纵轴并绘制第二个列表的数据
    ax2 = ax1.twinx()
    ax2.plot(x, CellProporty, 'bs--', label='Proporation')
    ax2.set_ylabel('Proportion in whole interactions')
    ax2.tick_params(axis='y')
    # 添加图例
    lines, labels = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines + lines2, labels + labels2)
    # 设置标题
    plt.title('Fig.1 Cells Chromatain Interner Count and Proporation')
    # 显示图形
    plt.savefig(save + '/Fig.1 Cells Chromatain Interner Count and Proporation.png',dpi=300)
    plt.close()


def DrawChromatainInternerCount(ChromIntraction = None, ChromLenList = None, save=''):

    ###Get Chrom Name
    ChromsName = set()
    for term in list(ChromIntraction[0].keys()):
        ChromsName.add(term.split("-")[0])
        ChromsName.add(term.split("-")[1])
    ChromsName = list(ChromsName)

    InteractionName = list(ChromIntraction[0].keys())

    InteractionArr = []
    for i in range(len(ChromIntraction)):
        InteractionArr.append([value for value in ChromIntraction[i].values()])
    InteractionArr = np.array(InteractionArr)

    ChromInner = []
    ChromOuter = [[] for _ in range(len(ChromsName))]

    for i in range(len(ChromsName)):
        for j in range(len(ChromsName)):

            if ChromsName[i] == ChromsName[j]:
                ChromInner.append(ChromsName[i]+"-"+ChromsName[j])
            else:
                ChromOuter[i].append(ChromsName[i]+"-"+ChromsName[j])
    ####inner
    ChromInnerpositions = []
    ChromInnerName = []
    for term in ChromInner:
        ChromInnerpositions.append(InteractionName.index(term))
        ChromInnerName.append(term.split("-")[0])
    ChromInnerInteraction = InteractionArr[:, ChromInnerpositions]

    ###Outner
    ChromOuterInteraction = []
    CellouterInteraction = []
    ChromOuterName = []
    for i in range(len(ChromsName)):
        curChromouterpositions = []
        ChromOuterName.append(ChromOuter[i][0].split("-")[0])

        for term in ChromOuter[i]:
            if term not in InteractionName:
                term = term.split("-")[1]+"-"+term.split("-")[0]
            curChromouterpositions.append(InteractionName.index(term))
        curChromOuterInteraction = InteractionArr[:, curChromouterpositions]
        CellouterInteraction.append(curChromOuterInteraction)
        curChromOuterInteraction = curChromOuterInteraction.flatten()
        ChromOuterInteraction.append(curChromOuterInteraction)


    #### Fig.2 Draw CellouterInteraction
    CellouterInteraction = np.hstack(CellouterInteraction)  #### CellouterInteraction
    CellouterInteraction = np.log10(CellouterInteraction + 1)
    # 计算每个细胞的平均值和标准差
    mean_values = np.mean(CellouterInteraction, axis=1)
    std_values = np.std(CellouterInteraction, axis=1)
    # sort
    mean_values = mean_values[np.argsort(mean_values)]
    std_values = std_values[np.argsort(mean_values)]
    # 创建一个图形对象和一个子图对象
    fig, ax = plt.subplots()
    # 绘制平均值线
    ax.errorbar(np.arange(len(mean_values)), mean_values, yerr=std_values, fmt='o', capsize=5)
    # 设置图形的标题和坐标轴标签
    ax.set_title('Fig.2 Cell Chromatain Outer Interactions')
    ax.set_xlabel('Cell number')
    ax.set_ylabel('1*10^N interactions')
    ax.set_ylim(bottom=0)
    # 显示图形
    plt.savefig(save + '/Fig.2 Cell Chromatain Outer Interactions.png',dpi=300)
    plt.close()


    #### Fig.3 Draw Chromatain Inner Interaction
    ChromInnerInteraction = np.log10(ChromInnerInteraction + 1)
    # 计算每个细胞的平均值和标准差
    mean_values = np.mean(ChromInnerInteraction, axis=0)
    std_values = np.std(ChromInnerInteraction, axis=0)
    # sort
    mean_values = mean_values[np.argsort(ChromsName)]
    std_values = std_values[np.argsort(ChromsName)]
    ChromsName_inner = [ChromsName[i] for i in np.argsort(ChromsName)] # 使用排序后的索引来获取对应的横坐标名字
    fig, ax = plt.subplots(figsize=(12, 12 / 1.618 ))
    # 绘制平均值线
    ax.errorbar(np.arange(len(mean_values)), mean_values, yerr=std_values, fmt='o', capsize=5)
    # 设置图形的标题和坐标轴标签
    x_ticks = np.arange(len(mean_values))
    x_tick_labels = ChromsName_inner
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(x_tick_labels)
    ax.set_title('Fig.3 Cell Chromatain Inner Interactions')
    ax.set_xlabel('Chromatain Names')
    ax.set_ylabel('1*10^N interactions')
    ax.set_ylim(bottom=0)
    # 显示图形
    plt.savefig(save + '/Fig.3 Cell Chromatain Inner Interactions.png',dpi=300)
    plt.close()

    #### Fig.4 Draw Chromatain Outer Interaction
    ChromOuterInteraction = np.array(ChromOuterInteraction)
    # 计算每个染色体的平均值和标准差
    mean_values = np.mean(ChromOuterInteraction, axis=1)
    std_values = np.std(ChromOuterInteraction, axis=1)
    # sort
    mean_values = mean_values[np.argsort(ChromsName)]
    std_values = std_values[np.argsort(ChromsName)]
    ChromsName_outer = [ChromsName[i] for i in np.argsort(ChromsName)] # 使用排序后的索引来获取对应的横坐标名字
    # 创建一个图形对象和一个子图对象
    fig, ax = plt.subplots(figsize=(12, 12 / 1.618 ))
    # 绘制平均值线
    ax.errorbar(np.arange(len(mean_values)), mean_values, yerr=std_values, fmt='o', capsize=5)
    # 设置图形的标题和坐标轴标签
    x_tick_labels = ChromsName_outer
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(x_tick_labels)
    ax.set_title('Fig.4 Chromatain Outer Interactions')
    ax.set_xlabel('Cell number')
    ax.set_ylabel('1*10^N interactions')
    ax.set_ylim(bottom=0)
    # 显示图形
    plt.savefig(save + '/Fig.4 Chromatain Outer Interactions.png',dpi=300)
    plt.close()

def DrawDistanceCellInternerCount(CellCounts=None, ChromLenList=None, save=None):
    # Fig.5 Cells Chromatain Interner distance interacions and Proporation
    ###Get Chrom Length
    values_list = [value for value in ChromLenList[0].values()]
    dataPoints = sum(np.multiply(values_list, values_list))
    ###Get propotion
    CellCounts = sorted(CellCounts)
    CellProporty = CellCounts / dataPoints
    ###Log10 count
    CellCounts = np.log10(np.array(CellCounts) + 1)
    ###Draw Plot
    x = list(range(1, len(CellCounts) + 1))
    # 创建图形和坐标轴对象
    fig, ax1 = plt.subplots()
    # 绘制第一个列表的数据
    ax1.plot(x, CellCounts, 'ro-', label='Count')
    ax1.set_xlabel('Cell Number')
    ax1.set_ylabel('1*10^n interactions')
    ax1.tick_params(axis='y')
    # 创建第二个纵轴并绘制第二个列表的数据
    ax2 = ax1.twinx()
    ax2.plot(x, CellProporty, 'bs--', label='Proporation')
    ax2.set_ylabel('Proportion in whole interactions')
    ax2.tick_params(axis='y')
    # 添加图例
    lines, labels = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines + lines2, labels + labels2)
    # 设置标题
    plt.title('Fig.5 Cells Chromatain Interner distance interacions and Proporation')
    # 显示图形
    plt.savefig(save + '/Fig.5 Cells Chromatain Interner distance interacions and Proporation.png', dpi=300)
    plt.close()


def DrawDistanceChromatainInternerCount(ChromatainCounts=None, save=None):
    # Fig.6 Chromatain distance Interactions
    InteractionArr = []
    for i in range(len(ChromatainCounts)):
        InteractionArr.append([value for value in ChromatainCounts[i].values()])
    InteractionArr = np.array(InteractionArr)
    ChromsName = list(ChromatainCounts[0].keys())

    ChromOuterInteraction = np.array(InteractionArr)
    # 计算每个染色体的平均值和标准差
    mean_values = np.mean(ChromOuterInteraction, axis=0)
    std_values = np.std(ChromOuterInteraction, axis=0)
    # sort
    mean_values = mean_values[np.argsort(ChromsName)]
    std_values = std_values[np.argsort(ChromsName)]
    ChromsNames = [ChromsName[i] for i in np.argsort(ChromsName)]  # 使用排序后的索引来获取对应的横坐标名字
    # 创建一个图形对象和一个子图对象
    fig, ax = plt.subplots(figsize=(12, 12 / 1.618))
    # 绘制平均值线
    ax.errorbar(np.arange(len(mean_values)), mean_values, yerr=std_values, fmt='o', capsize=5)
    # 设置图形的标题和坐标轴标签
    x_ticks = np.arange(len(mean_values))
    x_tick_labels = ChromsNames
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(x_tick_labels)
    ax.set_title('Fig.6 Chromatain distance Interactions')
    ax.set_xlabel('Cell number')
    ax.set_ylabel('1*10^N interactions')
    ax.set_ylim(bottom=0)
    # 显示图形
    plt.savefig(save + '/Fig.6 Chromatain distance Interactions.png', dpi=300)
    plt.close()



















