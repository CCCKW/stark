import os
import glob
import cooler
import numpy as np
from scipy.stats import norm
from .embedding import *
import pandas as pd

class MultiMcoolReader:
    def __init__(self,directory, resolution):
        self.directory = directory
        self.mcool_files = []
        self.coolers = []
        self.resolution = resolution
   

    def extract_mcool_files(self):
        self.mcool_files = []
        self.mcool_names = []
        self.GSE = []
        for root, dir_, files in os.walk(self.directory):
            for f in files:
                if f.endswith('.mcool'):
                    self.mcool_files.append(root + '/' + f)
                    self.mcool_names.append(f.split('.')[0])
                    expri = root.split('_' +f.split('.')[0]+ '_tmp')[0].split('/')[-1]
                    clss = root.split(expri+'_'+f.split('.')[0]+'_tmp')[0].split('/')[-2]
                    self.GSE.append(clss)


    def read_mcool_files(self):
        for mcool_file in self.mcool_files:
            try:
                cooler_obj = cooler.Cooler(mcool_file+"::/resolutions/"+str(self.resolution))
                self.coolers.append(cooler_obj)
            except Exception as e:
                print(f"Error reading file {mcool_file}: {str(e)}")

    def count_values(self):
        ###此函数为了计算每个single cell HIC matrix 有值的个数
        counts = []
        chromCounts = []
        for cooler_obj in self.coolers:
            count = 0
            currentChromDict = {}
            for chrom in cooler_obj.chromnames:
                matrix = cooler_obj.matrix(balance=False).fetch(chrom)
                chromCount = np.count_nonzero(matrix)
                count += chromCount
                currentChromDict[chrom] = chromCount
            counts.append(count)
            chromCounts.append(currentChromDict)

        return counts, chromCounts

    def get_coolers(self):
        return self.coolers


    def count_distance_interaction(self, RandomSize=3):
        ###此函数为了计算每个single cell HIC matrix 远端有值的个数，先得到矩阵，再生成一个相同大小、离对角线越近越可能是0、越远越可能是1的矩阵，之后再相乘
        counts = []
        chromCounts = []
        chromLens = []
        for cooler_obj in self.coolers:
            count = 0
            currentChromDict = {}
            currentChromLen = {}
            for chrom in cooler_obj.chromnames:
                matrix = cooler_obj.matrix(balance=False).fetch(chrom)
                matrixLen = np.shape(matrix)[0]
                distanceCount=0
                for i in range(RandomSize):
                    randomMat = self.generate_random_matrix(matrixLen)
                    matrix = matrix * randomMat
                    distanceCount += np.count_nonzero(matrix)
                count += distanceCount / RandomSize
                currentChromDict[chrom] = distanceCount / RandomSize
                currentChromLen[chrom] = matrixLen
            counts.append(count)
            chromCounts.append(currentChromDict)
            chromLens.append(currentChromLen)
        return counts, chromCounts, chromLens


    def generate_random_matrix(self, size):
        # Generate a matrix of random values between 0 and 1
        random_matrix = np.random.rand(size, size)

        # Calculate the distances from the main diagonal for all elements
        distances = np.abs(np.arange(size) - np.arange(size)[:, np.newaxis])

        # Calculate the probabilities for all elements based on distances
        probabilities = 1.0 - (distances / (size - 1))

        probabilities = probabilities * probabilities * probabilities
        # Generate the final matrix based on probabilities
        matrix = (random_matrix < probabilities).astype(int)

        return matrix

    def calculate_Chorm_interactions(self):
        ###计算染色体间的interaction
        interactions = []
        for cooler_obj in self.coolers:
            chromnames = cooler_obj.chromnames
            num_chroms = len(chromnames)
            chromnamesInteraction = {}
            for i in range(num_chroms):
                for j in range(i, num_chroms):
                    inter = cooler_obj.matrix(balance=False).fetch(chromnames[i], chromnames[j])
                    interaction = np.count_nonzero(inter)
                    chromnamesInteraction[chromnames[i]+"-"+chromnames[j]] = interaction
            interactions.append(chromnamesInteraction)
        return interactions

    def calculateAjcMatEmbadding(self, similarityMethod = 'gaussian_distance', embeddingMethod = 'tSNE', dim = 2,cache=True):
        CellSimilarity = []
        CellFeature = []
        fileCount = len(self.coolers)
        if not cache or not os.path.exists('./cache/cache.npz'):
            # 不存在cache，那还是要做一遍，但是存下来
            for i in range(fileCount):
                ####获取该细胞的HIC矩阵
                TargetMatrix_1 = self.coolers[i].matrix(balance=False)[:]

                ###获取距离对角线0-10的数据并降成一维
                diagonal_data = [np.diag(TargetMatrix_1, k) for k in range(0, 10)]

                ###获取降维后的矩阵
                SampledTargetMatrix_1 = self.scale_matrix(TargetMatrix_1, 200)

                ###整合整体特征与染色体内部特征
                CellFeature.append(np.concatenate([np.concatenate([arr.flatten() for arr in diagonal_data]), SampledTargetMatrix_1.flatten()], axis=0))

                ####计算细胞间的similarity
                CurCellSimilarity = []
                for j in range(fileCount):
                    TargetMatrix_2 = self.coolers[j].matrix(balance=False)[:]
                    SampledTargetMatrix_2 = self.scale_matrix(TargetMatrix_2, 200)
                    if similarityMethod == 'gaussian_distance':
                        CurCellSimilarity.append(self.gaussian_distance(SampledTargetMatrix_1, SampledTargetMatrix_2, 1))
                    if similarityMethod == 'EuclideanDistance':
                        CurCellSimilarity.append(np.linalg.norm(SampledTargetMatrix_1 - SampledTargetMatrix_2))
                    if similarityMethod == 'CosineSimilarity':
                        CurCellSimilarity.append(1 - self.cosine_similarity(SampledTargetMatrix_1, SampledTargetMatrix_2))
                ###记录每个细胞与所有细胞间的similarity
                CellSimilarity.append(CurCellSimilarity)
        else:
            ###存在cache，那就直接读取
            print('load cache')
            cache = np.load('./cache/cache.npz')
            CellFeature = cache['CellFeature']
            CellSimilarity = cache['CellSimilarity']

        
        
        # print(CellFeature)
        ##对矩阵正态分布处理
        CellSimilarity = np.array(CellSimilarity)
        smoothed_CellSimilarity = 1 - norm.cdf(CellSimilarity, np.mean(CellSimilarity), np.std(CellSimilarity))
        smoothed_CellSimilarity = (smoothed_CellSimilarity - np.min(smoothed_CellSimilarity)) / (np.max(smoothed_CellSimilarity) - np.min(smoothed_CellSimilarity))

        ###计算每个细胞的embadding
        CellFeature = np.array(CellFeature)
        self.CellFeature = CellFeature
        self.cellSimilarity = smoothed_CellSimilarity
        if cache and not os.path.exists('./cache/cache.npz'):
            if not os.path.exists('cache'):
                os.makedirs('cache')
            print('save cache')
            cache = {}
            cache['CellFeature'] = CellFeature
            cache['CellSimilarity'] = smoothed_CellSimilarity
            cache['mcool_names'] = self.mcool_names
            cache['GSE'] = self.GSE
            np.savez('./cache/cache.npz', **cache)

        if embeddingMethod == 'PCA':
            CellsEmbedding = PCAEmbedding(CellFeature, dim=dim)

        if embeddingMethod == 'MDS':
            CellsEmbedding = MDSEmbedding(CellFeature, dim=dim)

        if embeddingMethod == 'tSNE':
            CellsEmbedding = tSNEEmbedding(CellFeature, dim=dim)

        if embeddingMethod == 'UMAP':
            CellsEmbedding = UMAPEmbedding(CellFeature, dim=dim)

        if embeddingMethod == 'LLE':
            CellsEmbedding = LLEEmbedding(CellFeature, dim=dim)

        if embeddingMethod == 'Isomap':
            CellsEmbedding = IsomapEmbedding(CellFeature, dim=dim)

        return smoothed_CellSimilarity, CellsEmbedding

    

    def gaussian_distance(self, matrix1, matrix2, sigma):
        # 计算高斯核矩阵
        kernel1 = np.exp(-np.sum((matrix1[:, np.newaxis] - matrix1) ** 2, axis=2) / (2 * sigma ** 2))
        kernel2 = np.exp(-np.sum((matrix2[:, np.newaxis] - matrix2) ** 2, axis=2) / (2 * sigma ** 2))

        # 计算高斯距离
        distance = np.sum((kernel1 - kernel2) ** 2)

        return distance


    def cosine_similarity(self, matrix_A, matrix_B):
        ####计算两个矩阵的余弦相似性
        vector_A = matrix_A.ravel()  # 使用ravel()代替flatten()，无需创建副本
        vector_B = matrix_B.ravel()

        dot_product = np.dot(vector_A, vector_B)
        norm_A = np.linalg.norm(vector_A)
        norm_B = np.linalg.norm(vector_B)

        similarity = dot_product / (norm_A * norm_B)
        return similarity

    def scale_matrix(self, matrix, target_size):
        n = matrix.shape[0]  # 原始矩阵的大小
        scale_factor = n / target_size  # 计算缩放比例

        # 构建目标矩阵的行和列索引
        row_indices = np.arange(target_size) * scale_factor
        col_indices = np.arange(target_size) * scale_factor

        # 使用索引进行插值计算
        row_floor = row_indices.astype(int)
        row_ceil = np.ceil(row_indices).astype(int)
        col_floor = col_indices.astype(int)
        col_ceil = np.ceil(col_indices).astype(int)

        # 计算四个邻近像素值
        Q11 = matrix[row_floor[:, None], col_floor]
        Q12 = matrix[row_floor[:, None], col_ceil]
        Q21 = matrix[row_ceil[:, None], col_floor]
        Q22 = matrix[row_ceil[:, None], col_ceil]

        # 双线性插值计算目标矩阵
        target_matrix = (Q11 * (row_ceil[:, None] - row_indices[:, None]) * (col_ceil - col_indices) +
                         Q21 * (row_indices[:, None] - row_floor[:, None]) * (col_ceil - col_indices) +
                         Q12 * (row_ceil[:, None] - row_indices[:, None]) * (col_indices - col_floor) +
                         Q22 * (row_indices[:, None] - row_floor[:, None]) * (col_indices - col_floor))

        return target_matrix
