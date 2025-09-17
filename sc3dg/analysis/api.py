import pandas as pd
import os

import cooler
import numpy as np
from scipy.stats import norm
import sc3dg.cooltools.cooltools as cooltools
import warnings

warnings.filterwarnings("ignore")
from sc3dg.analysis.embedding import scale_matrix, gaussian_distance, cosine_similarity, PCAEmbedding, MDSEmbedding, tSNEEmbedding, \
    UMAPEmbedding, LLEEmbedding, IsomapEmbedding

from sc3dg.analysis.embedding import scale_matrix, gaussian_distance, cosine_similarity, PCAEmbedding, MDSEmbedding, tSNEEmbedding, \
    UMAPEmbedding, LLEEmbedding, IsomapEmbedding
from sc3dg.utils.gini import calculate_giniqc

import bioframe
import sys

from sc3dg.cooltools.cooltools import insulation
import random
import tqdm

import matplotlib.pyplot as plt
import seaborn as sns

from joblib import Parallel, delayed
hg38_cens = pd.DataFrame({'chrom': ['chr1', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17',
                                    'chr18', 'chr19', 'chr2', 'chr20', 'chr21', 'chr22', 'chr3', 'chr4', 'chr5', 'chr6',
                                    'chr7', 'chr8', 'chr9', 'chrX', 'chrY'],
                          'start': [122026459, 39686682, 51078348, 34769407, 16000000, 16000000, 17083673, 36311158,
                                    22813679, 15460899, 24498980, 92188145, 26436232, 10864560, 12954788, 90772458,
                                    49712061, 46485900, 58553888, 58169653, 44033744, 43389635, 58605579, 10316944],
                          'end': [124932724, 41593521, 54425074, 37185252, 18051248, 18173523, 19725254, 38265669,
                                  26616164, 20861206, 27190874, 94090557, 30038348, 12915808, 15054318, 93655574,
                                  51743951, 50059807, 59829934, 61528020, 45877265, 45518558, 62412542, 10544039],
                          'mid': [123479591, 40640101, 52751711, 35977329, 17025624, 17086761, 18404463, 37288413,
                                  24714921, 18161052, 25844927, 93139351, 28237290, 11890184, 14004553, 92214016,
                                  50728006, 48272853, 59191911, 59848836, 44955504, 44454096, 60509060, 10430491]})

mm10_cens = pd.DataFrame(columns=['chrom', 'start', 'end', 'mid'])
mm10_cens['start'] = mm10_cens['start'].astype('int64')
mm10_cens['end'] = mm10_cens['end'].astype('int64')
mm10_cens['mid'] = mm10_cens['mid'].astype('int64')

mm10_chromsize = pd.DataFrame({'chrom': ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
                                         'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17',
                                         'chr18', 'chr19', 'chrX', 'chrY', 'chrM'],
                               'length': [195471971, 182113224, 160039680, 156508116, 151834684, 149736546, 145441459,
                                          129401213, 124595110, 130694993, 122082543, 120129022, 120421639, 124902244,
                                          104043685, 98207768, 94987271, 90702639, 61431566, 171031299, 91744698,
                                          16299]})
hg38_chromsize = pd.DataFrame({'chrom': ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
                                         'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17',
                                         'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM'],
                               'length': [248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973,
                                          145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718,
                                          101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983,
                                          50818468, 156040895, 57227415, 16569]})
data_info = {'hg38': [hg38_chromsize, hg38_cens], 'mm10': [mm10_chromsize, mm10_cens]}


class MMCooler:
    def __init__(self,
                 directory,  # 目录
                 resolution,  # 分辨率
                 describe=None,  # 技术名称
                 merge=False,  # 是否将读取的cell全部merge
                 genome: str = 'mm10',
                 balance: bool = True,
                 read_cache=True,
                 n_jobs=6,
                 ):
        '''
            directory:
            resolution:如果int,那没什么，就是直接读取；如果是list，list第一个必须要有，切是最小的，后面的分辨率会按照这个去合并或者读取
            describe:
            merge:
            genome:
            balance
        '''
        self.directory = directory
        self.n_jobs = n_jobs
        self.mcool_files = []
        self.mcool_names = []
        self.cache_file = []
        self.clr_obj = None

        self.__merge = merge
        self.genome = genome
        self.__balance = balance
       
        self.resolution = resolution if isinstance(resolution, list) else [resolution]
        self.resolution = [int(x) for x in resolution]
        self.describe = describe
        self.chromosome = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6',
                           'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12',
                           'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18',
                           'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM']
        self.__chrom_terr = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6',
                             'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12',
                             'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18',
                             'chr19']
        self.read_cache = read_cache
        self.__pool = False

        if not os.path.exists('./cache'):
            os.makedirs('./cache')
        self.extract_mcool_files()
        self.read_mcool_files()
        if not os.path.exists('./file'):
            os.makedirs('./file')

    def stat_merge(self):
        return self.__merge

    def stat_balance(self):
        return self.__balance

    def stat_resolution(self):
        return isinstance(self.resolution, list)

    def shape(self):
        shape = np.array(self.clr_obj)
        return shape

    def stat_file_name(self):
        return self.mcool_names

    def read_zoomify(self, mcool_file, cache_name):
        raw_res = self.resolution[0] 
        if not os.path.exists(cache_name):
            cooler.zoomify_cooler(base_uris=mcool_file + "::resolutions/" + str(raw_res), outfile=cache_name,
                                    resolutions=self.resolution,
                                    chunksize=10000, nproc=16)

    # 读取目录下所有的mcool文件
    def extract_mcool_files(self):
        if self.directory.endswith('.mcool'):
            self.mcool_files.append(self.directory)
            self.mcool_names.append(self.directory.split('/')[-1].split('.')[0])

        for root, dir_, files in os.walk(self.directory):
            for f in files:
                if f.endswith('.mcool'):
                    self.mcool_files.append(root + '/' + f)
                    self.mcool_names.append(f.split('.')[0])

    # 加载所有mcool文件到cooler.Cooler对象,如果merge不是none，那就merge到这个
    def read_mcool_files(self):
        '''
            self.coolers: list of cooler objects
        '''
        # 如果合并
        if self.__merge:
            print('run merge mode')
            # 多个分辨率
            self.clr_obj = []
            raw_res = self.resolution[0]  # 读取有这个分辨率,并以此分辨率座位基础

            # 先把每个样本bin化，然后再读取
            print('binning the raw mcool files')
            for mcool_file, mcool_name in zip(self.mcool_files, self.mcool_names):
                # 如果美有cache，那就创建，读取
                # 创建对应样本cache的名字
                cache_name = './cache/%s_multi_res_%s' % (self.describe, mcool_name)
                for res in self.resolution:
                    cache_name += '-'
                    cache_name += str(res)
                cache_name += '.mcool'

                # 如果没有这个文件，那就生成
                if not os.path.exists(cache_name):

                    cooler.zoomify_cooler(base_uris=mcool_file + "::resolutions/" + str(raw_res), outfile=cache_name,
                                          resolutions=self.resolution,
                                          chunksize=10000, nproc=32)

            # 都存到cache了，然后按照分辨率读取,并 merge
            print('reading the binned mcool files')
            for res in self.resolution:
                tmp = []
                for mcool_file, mcool_name in zip(self.mcool_files, self.mcool_names):
                    cache_name = './cache/%s_multi_res_%s' % (self.describe, mcool_name)
                    for rr in self.resolution:
                        cache_name += '-'
                        cache_name += str(rr)
                    cache_name += '.mcool'
                    tmp.append(cache_name + '::resolutions/' + str(res))

                print('merging the binned mcool files/resolution=%s' % res)
                if not os.path.exists('./cache/%s_merge_all_res-%s' % (self.describe, res)):
                    cooler.merge_coolers(input_uris=tmp,
                                         output_uri='./cache/%s_merge_all_res-%s' % (self.describe, res),
                                         mergebuf=1000000)
                self.clr_obj.append(cooler.Cooler('./cache/%s_merge_all_res-%s' % (self.describe, res)))

        # 如果不合并
        else:
            print('run nomerge mode')
            self.clr_obj = []
            # 多个分辨率

            # 先把每个样本bin化，然后再读取
            raw_res = self.resolution[0]  # 读取有这个分辨率,并以此分辨率座位基础
            print('binning the raw mcool files')
            task = []
            for mcool_file, mcool_name in zip(self.mcool_files, self.mcool_names):
                cache_name = './cache/%s_multi_nomerge_%s_res' % (self.describe, mcool_name)
                for res in self.resolution:
                    cache_name += '-'
                    cache_name += str(res)
                cache_name += '.mcool'
                task.append([ mcool_file,  cache_name])

            Parallel(n_jobs=self.n_jobs)(delayed(self.read_zoomify)(mcool_file, cache_name) for mcool_file, cache_name in task)
            
            # 都存到cache了，然后按照分辨率读取
            # 外层是分辨率，里层是样本
            print('reading the binned mcool files')
            for res in self.resolution:
                print(res)
                tmp = []
                for mcool_file, mcool_name in zip(self.mcool_files, self.mcool_names):
                    cache_name = './cache/%s_multi_nomerge_%s_res' % (self.describe, mcool_name)
                    for r in self.resolution:
                        cache_name += '-'
                        cache_name += str(r)
                    cache_name += '.mcool'
                    tmp.append(cooler.Cooler(cache_name + '::resolutions/' + str(res)))
                    if cache_name not in self.cache_file:
                        self.cache_file.append(cache_name)
                self.clr_obj.append(tmp)

    # 获取单个细胞所有contact数量
    def get_total_contact(self, res=None):
        '''
            res: resolution
        '''
        if res is None:
            clr = self.clr_obj
        else:
            clr = [self.clr_obj[self.resolution.index(res)]]

        res = []
        for cr in clr[0]:
            res.append(cr.info['sum'])
        df = pd.DataFrame(res, index=self.mcool_names, columns=['total_contact'])
        df.to_csv('./file/{}_total_contact.csv'.format(self.describe), sep='\t')
        return res

    # 输入一个cooler，返回一组结果
    def __get_contact_vs_distance(self, clr, res, arms):
        resolution = res
        arms = arms[arms.chrom.isin(clr.chromnames)].reset_index(drop=True)
        if not cooltools.lib.is_cooler_balanced(clr):
            cooler.balance_cooler(clr, store=True)
        cvd_smooth_agg = cooltools.expected_cis(
            clr=clr,
            view_df=arms,
            smooth=True,
            aggregate_smoothed=True,
            nproc=16
        )

        cvd_smooth_agg['s_bp'] = cvd_smooth_agg['dist'] * resolution
        cvd_smooth_agg['balanced.avg.smoothed.agg'].loc[cvd_smooth_agg['dist'] < 2] = np.nan
        cvd_smooth_agg = cvd_smooth_agg.groupby('s_bp').sum().reset_index()
        cvd_smooth_agg = cvd_smooth_agg.set_index('s_bp')

        return cvd_smooth_agg.loc[:, ['balanced.avg.smoothed.agg']]

    # 获得绘制P（s）图像的数据
    def get_contact_vs_distance(self, res=None, version='hg38'):
        '''
            dismiss:舍弃对角线开始的第几个细胞
            mode:
                cis:   interaction decay with distance
                 or trans:teraction frequencies for inter-chromosomal blocks
            intra_only:  Return expected only for symmetric intra-regions defined by view_df,
            i.e. chromosomes, chromosomal-arms, intra-domains, etc.
            When False returns expected both for symmetric intra-regions and
            assymetric inter-regions.
        '''

        if res is None:
            sys.exit('res is None')

        chromsizes, cens = data_info[version]
        arms = bioframe.make_chromarms(chromsizes, cens)

        clr = self.clr_obj[self.resolution.index(res)]
        result = []
        for name, val in zip(self.mcool_names, clr):
            tmp = self.__get_contact_vs_distance(val, res, arms)
            tmp.columns = [name]
            result.append(tmp)
        result = pd.concat(result, axis=1)

        return result

    # 获得intra-chrom contact个数
    def get_intra_contact_count(self, split_chorm: bool = True):
        '''
            merge + split: 返回染色体dict
            merge + nosolt:返回一个值
            nomerge + nosplit:返回list of count
            nomerge + split:返回 list of dict
        '''

        counts = []

        for val in self.clr_obj[-1]:
            counts.append(self.__get_intra_contact_count(val, split_chorm))

        df = pd.DataFrame(counts)
        df.index = self.mcool_names

        df.to_csv('./file/{}_splitChrom:{}_intra_contact_count.csv'.format(self.describe, split_chorm), sep='\t')

        return counts

    # 一个cool，一个dict
    def __get_intra_contact_count(self, clr, split_chorm: bool = True):

        bins = clr.bins()[:]
        pixels = clr.pixels()[:]
        pixels['chr1'] = list(bins.iloc[list(pixels['bin1_id']), 0])
        pixels['chr2'] = list(bins.iloc[list(pixels['bin2_id']), 0])
        pixels = pixels[pixels['chr1'] == pixels['chr2']]
        if split_chorm:
            # 如果merge，那就返回一个dict
            counts = {}
            for chrom in clr.chromnames:
                counts[chrom] = pixels[pixels['chr1'] == chrom].shape[0]
        # 不分染色体，那就只返回一个数值
        else:
            counts = {'total': 0}
            for chrom in clr.chromnames:
                counts['total'] += pixels[pixels['chr1'] == chrom].shape[0]

        return counts

    # 获得contact agg,只算上三角+对角线，求和
    def get_intra_contact_agg(self, split_chorm: bool = True, res=None):
        '''
            merge + split: 返回染色体dict
            merge + nosolt:返回一个值
            nomerge + nosplit:返回list of count
            nomerge + split:返回 list of dict
        '''
        counts = []

        for val in self.clr_obj[-1]:
            counts.append(self.__get_intra_contact_agg(val, split_chorm))

        df = pd.DataFrame(counts)
        df.index = self.mcool_names

        df.to_csv('./file/{}_splitChrom:{}_intra_contact_agg.csv'.format(self.describe, split_chorm), sep='\t')
        return counts

    # 一个cooelr，一个result
    def __get_intra_contact_agg(self, clr, split_chorm: bool = True):
        # 分染色体计数
        bins = clr.bins()[:]
        pixels = clr.pixels()[:]
        pixels['chr1'] = list(bins.iloc[list(pixels['bin1_id']), 0])
        pixels['chr2'] = list(bins.iloc[list(pixels['bin2_id']), 0])
        pixels['start'] = list(bins.iloc[list(pixels['bin1_id']), 1])
        pixels['end'] = list(bins.iloc[list(pixels['bin1_id']), 2])
        pixels = pixels[pixels['chr1'] == pixels['chr2']]
        if split_chorm:
            # 如果merge，那就返回一个dict

            counts = {}
            for chrom in self.chromosome:
                try:
                    tmp = pixels[(pixels['chr1'] == chrom)]
                except:
                    continue

                counts[chrom] = sum(tmp['count'])
        # 不分染色体，那就只返回一个数值
        else:

            counts = {'total': sum(pixels[pixels['chr1'] == pixels['chr2']]['count'])}

        return counts

    # 获得inter-chrom contact个数
    def get_inter_contact_count(self):
        counts = []

        for val in self.clr_obj[-1]:
            counts.append(self.__get_inter_contact_count(val))
        df = pd.DataFrame(counts)
        df.index = self.mcool_names

        df.to_csv('./file/{}_inter_contact_count.csv'.format(self.describe), sep='\t')
        return counts

    def __get_inter_contact_count(self, clr):
        # 分染色体计数
        bins = clr.bins()[:]
        pixels = clr.pixels()[:]
        pixels['chr1'] = list(bins.iloc[list(pixels['bin1_id']), 0])
        pixels['chr2'] = list(bins.iloc[list(pixels['bin2_id']), 0])
        pixels['start'] = list(bins.iloc[list(pixels['bin1_id']), 1])
        pixels['end'] = list(bins.iloc[list(pixels['bin1_id']), 2])
        pixels = pixels[pixels['chr1'] != pixels['chr2']]
        counts = {'total': pixels.shape[0]}

        return counts

    # 获得inter-chrom contact总数
    def get_inter_contact_agg(self, split_chrom: bool = True, res=None):
        counts = []

        for val in self.clr_obj[-1]:
            counts.append(self.__get_inter_contact_agg(val))
        df = pd.DataFrame(counts)
        df.index = self.mcool_names

        df.to_csv('./file/{}_inter_contact_agg.csv'.format(self.describe), sep='\t')
        return counts

    def __get_inter_contact_agg(self, clr, split_chorm: bool = True):
        # 分染色体计数
        bins = clr.bins()[:]
        pixels = clr.pixels()[:]
        pixels['chr1'] = list(bins.iloc[list(pixels['bin1_id']), 0])
        pixels['chr2'] = list(bins.iloc[list(pixels['bin2_id']), 0])
        pixels['start'] = list(bins.iloc[list(pixels['bin1_id']), 1])
        pixels['end'] = list(bins.iloc[list(pixels['bin1_id']), 2])
        pixels = pixels[pixels['chr1'] != pixels['chr2']]
        counts = {'total': sum(pixels['count'])}

        return counts

    # 获得clr_obj对象中的embeeding,仅在merge=False时候使用
    def get_clustering_embedding(self, res=None, cache=True,
                                 similarityMethod='gaussian_distance',
                                 embeddingMethod='UMAP',
                                 dim=2):
        if self.__merge:
            print('merged,no embedding')
            return

        if res is None:
            print('please input resolution')
            return

        CellSimilarity = []
        CellFeature = []

        clr = self.clr_obj[self.resolution.index(res)]

        fileCount = len(clr)

        if not cache or not os.path.exists('./cache/cache_%s_%s.npz' % (self.describe, res)):
            # 不存在cache，那还是要做一遍，但是存下来

            for cr in clr:
                ####获取该细胞的HIC矩阵
                TargetMatrix_1 = cr.matrix(balance=False)[:]

                ###获取距离对角线0-10的数据并降成一维
                diagonal_data = [np.diag(TargetMatrix_1, k) for k in range(0, 10)]

                ###获取降维后的矩阵
                SampledTargetMatrix_1 = scale_matrix(TargetMatrix_1, 200)

                ###整合整体特征与染色体内部特征
                CellFeature.append(np.concatenate(
                    [np.concatenate([arr.flatten() for arr in diagonal_data]), SampledTargetMatrix_1.flatten()],
                    axis=0))

                ####计算细胞间的similarity
                CurCellSimilarity = []
                for j in range(fileCount):
                    TargetMatrix_2 = clr[j].matrix(balance=False)[:]
                    SampledTargetMatrix_2 = scale_matrix(TargetMatrix_2, 200)
                    if similarityMethod == 'gaussian_distance':
                        CurCellSimilarity.append(gaussian_distance(SampledTargetMatrix_1, SampledTargetMatrix_2, 1))
                    if similarityMethod == 'EuclideanDistance':
                        CurCellSimilarity.append(np.linalg.norm(SampledTargetMatrix_1 - SampledTargetMatrix_2))
                    if similarityMethod == 'CosineSimilarity':
                        CurCellSimilarity.append(1 - cosine_similarity(SampledTargetMatrix_1, SampledTargetMatrix_2))
                ###记录每个细胞与所有细胞间的similarity
                CellSimilarity.append(CurCellSimilarity)
        else:
            ###存在cache，那就直接读取
            print('load cache')
            cache = np.load('./cache/cache_%s_%s.npz' % (self.describe, res))
            CellFeature = cache['CellFeature']
            CellSimilarity = cache['CellSimilarity']

        # print(CellFeature)
        ##对矩阵正态分布处理
        CellSimilarity = np.array(CellSimilarity)
        smoothed_CellSimilarity = 1 - norm.cdf(CellSimilarity, np.mean(CellSimilarity), np.std(CellSimilarity))
        smoothed_CellSimilarity = (smoothed_CellSimilarity - np.min(smoothed_CellSimilarity)) / (
                    np.max(smoothed_CellSimilarity) - np.min(smoothed_CellSimilarity))

        ###计算每个细胞的embadding
        CellFeature = np.array(CellFeature)
        self.CellFeature = CellFeature
        self.cellSimilarity = smoothed_CellSimilarity

        if cache and not os.path.exists('./cache/cache_%s_%s.npz' % (self.describe, res)):
            if not os.path.exists('cache'):
                os.makedirs('cache')
            print('save cache')
            cache = {}
            cache['CellFeature'] = CellFeature
            cache['CellSimilarity'] = smoothed_CellSimilarity
            cache['mcool_names'] = self.mcool_names

            np.savez('./cache/cache_%s_%s.npz' % (self.describe, res), **cache)

        print(CellFeature)
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
        else:
            print('please input right embedding method')
            return smoothed_CellSimilarity

        return smoothed_CellSimilarity, CellsEmbedding

    # 获得每个cell的对角线元素flateen的feature
    def get_clustering_diagnal(self, dim=10, res=None):

        shape = np.array(self.clr_obj).shape

        if len(shape) == 1:
            CellFeature = []
            for clr in self.clr_obj:
                ####获取该细胞的HIC矩阵
                TargetMatrix_1 = clr.matrix(balance=False)[:]

                ###获取距离对角线0-10的数据并降成一维
                diagonal_data = [np.diag(TargetMatrix_1, k) for k in range(0, dim)]

                ###获取降维后的矩阵
                SampledTargetMatrix_1 = scale_matrix(TargetMatrix_1, 200)

                ###整合整体特征与染色体内部特征
                result = np.concatenate(
                    [np.concatenate([arr.flatten() for arr in diagonal_data]), SampledTargetMatrix_1.flatten()], axis=0)
                CellFeature.append(result)
        else:
            CellFeature = []
            for lt in self.clr_obj:
                tmp_feature = []
                for clr in lt:
                    ####获取该细胞的HIC矩阵
                    TargetMatrix_1 = clr.matrix(balance=True)[:]

                    ###获取距离对角线0-10的数据并降成一维
                    diagonal_data = [np.diag(TargetMatrix_1, k) for k in range(0, dim)]

                    ###获取降维后的矩阵
                    SampledTargetMatrix_1 = scale_matrix(TargetMatrix_1, 200)

                    ###整合整体特征与染色体内部特征
                    result = np.concatenate(
                        [np.concatenate([arr.flatten() for arr in diagonal_data]), SampledTargetMatrix_1.flatten()],
                        axis=0)
                    tmp_feature.append(result)
                CellFeature.append(tmp_feature)
        return CellFeature

    # 获得每个mcool文件的contact per bin
    def get_avg_contact_per_bin(self, ignore_diag=False, res=None):
        if res is None:
            clr = self.clr_obj
        else:
            clr = self.clr_obj[self.resolution.index(res)]

        shape = np.array(clr).shape

        if len(shape) == 1:
            avgs = []
            for cr in clr:
                avgs.append(self.__get_avg_contact_per_bin(cr, ignore_diag))
        else:
            avgs = []
            for val in clr:
                tmp_avg = []
                for cr in val:
                    tmp_avg.append(self.__get_avg_contact_per_bin(cr, ignore_diag))
                avgs.append(tmp_avg)

        return avgs

    def __get_avg_contact_per_bin(self, clr, ignore_diag=False):
        top = 0
        down = 0
        for chrom in self.chromosome:

            try:
                mat = clr.matrix(balance=False).fetch(chrom)
            except:
                continue
            if ignore_diag:
                top += np.triu(mat, 1).sum()
            else:
                top += np.triu(mat).sum()

            down += mat.shape[0]

        return top / down

    # 获得不同分辨的覆盖率
    def get_coverage(self, res=None):

        if res is None:
            clr = self.clr_obj
        else:
            clr = [self.clr_obj[self.resolution.index(res)]]

        coverages = []
        for val in clr:
            tmp_coverage = []
            for cr in val:
                tmp_coverage.append(self.__get_coverage(cr))
            coverages.append(tmp_coverage)

        if res is None:
            for res in self.resolution:
                df = pd.DataFrame(coverages[self.resolution.index(res)])
                df.index = self.mcool_names
                df.to_csv('./file/{}_{}_coverages.csv'.format(self.describe, res), sep='\t')

        else:
            df = pd.DataFrame(coverages[0])
            df.index = self.mcool_names
            df.to_csv('./file/{}_{}_coverages.csv'.format(self.describe, res), sep='\t')

        return coverages

    def __get_coverage(self, clr):

        bins = clr.bins()[:]
        pixels = clr.pixels()[:]
        cover = set(set(np.unique(pixels['bin1_id'])).union(set(np.unique(pixels['bin2_id']))))
        freq = len(cover) / bins.shape[0]
        return freq

    # 获得累计覆盖率
    def get_acum_coverage(self, res=None):
        if self.__merge:
            sys.exit('merge is True, can not get acum_coverage')

        if res is None:
            clr = self.clr_obj
        else:
            clr = [self.clr_obj[self.resolution.index(res)]]

        coverages = []
        for val in clr:
            tmp_coverage = []
            for cr in val:
                tmp_coverage.append(self.__get_acum_coverage(cr))
            coverages.append(tmp_coverage)

        return coverages

    def __get_acum_coverage(self, clr):
        freq = []
        bins = clr[0].bins()[:]
        pixels = clr[0].pixels()[:]
        cover = set(set(np.unique(pixels['bin1_id'])).union(set(np.unique(pixels['bin2_id']))))
        freq.append(len(cover) / bins.shape[0])
        for cr in clr[1:]:
            bins = cr.bins()[:]
            pixels = cr.pixels()[:]
            cover = cover.union(set(np.unique(pixels['bin1_id'])).union(set(np.unique(pixels['bin2_id']))))
            freq.append(len(cover) / bins.shape[0])
        freq = sorted(freq)
        return freq

    def get_nDS_chrom_terr(self, res=None, epoch=100):
        if self.__merge:
            sys.exit('merge is True, can not get nDS_chrom_terr_per_cell')

        if res is None:
            clr = self.clr_obj
        else:
            clr = [self.clr_obj[self.resolution.index(res)]]

        self.chrom_pair = []
        for cr in self.__chrom_terr:
            for cr2 in self.__chrom_terr[self.__chrom_terr.index(cr) + 1:]:
                self.chrom_pair.append([cr, cr2])

        result = []
        for i, val in enumerate(clr):
            # print(i)
            result.append(self.__get_nDS_chrom_terr(val, i, epoch))

        return result

        # 获得单个细胞的每个pair的nDS和result

    def __get_nDS_chrom_terr(self, clr, idx, epoch):
        print('cal DS')
        DS = []
        for cell in clr:
            bins = cell.bins()[:]
            pair_DS = []
            pixels = cell.pixels()[:]
            # print(bins)
            pixels['chr1'] = list(bins.iloc[list(pixels['bin1_id']), 0])
            pixels['chr2'] = list(bins.iloc[list(pixels['bin2_id']), 0])
            for pair in self.chrom_pair:
                chr1 = list(bins[bins['chrom'] == pair[0]].index)
                chr2 = list(bins[bins['chrom'] == pair[1]].index)
                pair_DS.append(self.__terri_detection_score(pixels, chr1, chr2))
            DS.append(pair_DS)
        df = pd.DataFrame(DS)
        df.columns = [x[0] + '_' + x[1] for x in self.chrom_pair]

        df.to_csv('./file/%s_DS_%s.csv' % (self.describe, str(self.resolution[idx])), sep='\t')
        print('cal random DS')
        # random ds
        random_DS = []
        for cell in clr:

            pair_DS = [0] * 171
            for i in tqdm.tqdm(range(epoch)):
                pixels = cell.pixels()[:]
                for pair in self.chrom_pair:
                    chr1 = list(bins[bins['chrom'] == pair[0]].index)
                    chr2 = list(bins[bins['chrom'] == pair[1]].index)

                    tmp = pixels.sample(len(chr1) + len(chr2), replace=True)
                    pair_DS[self.chrom_pair.index(pair)] += self.__terri_detection_score(tmp, chr1, chr2)
            random_DS.append([x / epoch for x in pair_DS])
        df2 = pd.DataFrame(random_DS)
        df2.columns = [x[0] + '_' + x[1] for x in self.chrom_pair]

        df2.to_csv('./file/%s_terr_rDS_%s.csv' % (self.describe, str(self.resolution[idx])), sep='\t', index=None)

        nDS = df - df2

        nDS.to_csv('./file/%s_terr_nDS_%s.csv' % (self.describe, str(self.resolution[idx])), sep='\t', index=None)
        return nDS

    def __terri_detection_score(self, pixels: pd.DataFrame, chr1: list, chr2: list):
        AA = pixels[pixels['bin1_id'].isin(chr1) & pixels['bin2_id'].isin(chr1)].shape[0]
        BB = pixels[pixels['bin1_id'].isin(chr2) & pixels['bin2_id'].isin(chr2)].shape[0]
        AB = pixels[pixels['bin1_id'].isin(chr1) & pixels['bin2_id'].isin(chr2)].shape[0]
        BA = pixels[pixels['bin1_id'].isin(chr2) & pixels['bin2_id'].isin(chr1)].shape[0]
        Eintra = len(chr1) * len(chr1) + len(chr2) * len(chr2)
        Einter = len(chr1) * len(chr2) * 2
        # print(AA,BB,AB,BA,Eintra,Einter)
        DS = (AA + BB) / Eintra - (AB + BA) / Einter
        return DS

    # 计算每个bin的A/B compartment
    def get_AB_compartment(self, genome='', res=None, epoch=100):
        if self.__merge:
            sys.exit('merge is True, can not get AB_compartmen')
        if res is None:
            sys.exit('res is None, can not get AB_compartmen')
        if res not in self.resolution:
            sys.exit('res is not in resolution, can not get AB_compartmen')

        clr = self.clr_obj[self.resolution.index(res)]

        print('merging cooler')
        merge = []
        for val in self.cache_file:
            name = val + '::resolutions/' + str(res)
            merge.append(name)
        merge_name = './cache/merge_' + str(res) + '_' + self.describe + '.mcool'
        if not os.path.exists(merge_name):
            cooler.merge_coolers(input_uris=merge, output_uri=merge_name, mergebuf=1000000)
        tmp_merge = cooler.Cooler(merge_name)
        cooler.balance_cooler(tmp_merge, store=True)
        print('done')
        genome = bioframe.load_fasta(genome)
        bins = tmp_merge.bins()[:]
        gc_cov = bioframe.frac_gc(bins[['chrom', 'start', 'end']], genome)
        view_df = pd.DataFrame({'chrom': tmp_merge.chromnames,
                                'start': 0,
                                'end': tmp_merge.chromsizes.values,
                                'name': tmp_merge.chromnames}
                               )
        cis_eigs = cooltools.eigs_cis(
            tmp_merge,
            gc_cov,
            view_df=view_df,
            n_eigs=3,
        )

        # cis_eigs[0] returns eigenvalues, here we focus on eigenvectors
        eigenvector_track = cis_eigs[1][['chrom', 'start', 'end', 'E1']]
        store_name = 'E1'
        with tmp_merge.open("r+") as grp:
            if store_name in grp["bins"]:
                del grp["bins"][store_name]
            h5opts = {"compression": "gzip", "compression_opts": 6}
            grp["bins"].create_dataset(store_name, data=eigenvector_track['E1'], **h5opts)

        data = tmp_merge.bins()[:]

        data['E1'] = np.nan_to_num(data['E1'])
        data['sign'] = '$'
        data['sign'][data['E1'] > 0] = 'A'
        data['sign'][data['E1'] < 0] = 'B'
        a_list = []
        b_list = []
        for c in set(data['chrom']):
            seq = list(data[data['chrom'] == c]['sign'])
            base = data[data['chrom'] == c].index[0]

            seq = ''.join(seq)
            #     print(seq)
            new_seq = ''
            idx = []
            tmp_idx = []
            for i in range(len(seq)):
                sym = seq[i]
                if new_seq == '':
                    new_seq += sym
                    tmp_idx.append(i + base)
                else:
                    if sym == seq[i - 1]:
                        tmp_idx.append(i + base)
                    else:
                        idx.append(tmp_idx)
                        tmp_idx = [i + base]
                        new_seq += sym
            for i in range(len(new_seq) - 3):
                if new_seq[i:i + 3] == 'ABA':
                    tmp_A = idx[i] + idx[i + 2]
                    tmp_B = idx[i + 1]
                    a_list.append(tmp_A)
                    b_list.append(tmp_B)

                elif new_seq[i:i + 3] == 'BAB':
                    tmp_A = idx[i + 1]
                    tmp_B = idx[i] + idx[i + 2]
                    a_list.append(tmp_A)
                    b_list.append(tmp_B)
        result = []
        for i in tqdm.tqdm(range(len(clr))):
            cr = clr[i]
            bins = cr.bins()[:]
            pixels = cr.pixels()[:]
            pixels['chr1'] = list(bins.iloc[list(pixels['bin1_id']), 0])
            pixels['chr2'] = list(bins.iloc[list(pixels['bin2_id']), 0])
            pixels['start'] = list(bins.iloc[list(pixels['bin1_id']), 1])
            pixels['end'] = list(bins.iloc[list(pixels['bin1_id']), 2])
            pixels = pixels[pixels['chr1'] == pixels['chr2']]
            DS = []
            for i, j in zip(a_list, b_list):
                DS.append(self.__terri_detection_score(pixels, i, j))
            rDS = self.__return_AB_rds(a_list, b_list, pixels, bins, epoch)
            result.append((DS - rDS).mean())

        return result

    def __return_AB_rds(self, a_list, b_list, pixels, bins, epoch):

        pixels['chr1'] = list(bins.iloc[list(pixels['bin1_id']), 0])
        pixels['chr2'] = list(bins.iloc[list(pixels['bin2_id']), 0])
        pixels['start'] = list(bins.iloc[list(pixels['bin1_id']), 1])
        pixels['end'] = list(bins.iloc[list(pixels['bin1_id']), 2])
        pixels = pixels[pixels['chr1'] == pixels['chr2']]
        bins['index'] = list(bins.index)
        bins = bins.groupby('chrom').describe()['index'].loc[:, ['min', 'max']].astype(int)
        res = []
        for _ in range(epoch):
            for i, j in zip(a_list, b_list):
                idx = np.random.randint(bins['min'], bins['max'])
                idx = random.sample(list(idx), 1)[0]
                tmp = list(sorted(i + j))

                if tmp.index(j[0]) != 0:
                    ch1 = list(range(idx, idx + tmp.index(j[0]))) + list(
                        range(idx + tmp.index(j[-1]) + 1, idx + len(tmp)))
                    ch2 = list(range(idx + tmp.index(j[0]), idx + tmp.index(j[-1]) + 1))
                else:
                    ch1 = list(range(idx, idx + tmp.index(j[-1]) + 1))
                    ch2 = list(range(idx + tmp.index(j[-1]) + 1, idx + len(tmp)))
                try:
                    rds = self.__terri_detection_score(pixels, ch1, ch2)
                except:
                    continue
                res.append(rds)
        return np.array(res).mean()

    def get_TAD_nDS(self, res=500000, epoch=500, nbins=10, ):
        chrom = ['chr' + str(i) for i in range(20)]
        if res not in self.resolution:
            sys.exit('resolution not in %s' % str(self.resolution))

        merge = []
        for val in self.cache_file:
            name = val + '::resolutions/' + str(res)
            merge.append(name)
        merge_name = './cache/merge_' + str(res) + '_' + self.describe + '.mcool'
        if not os.path.exists(merge_name):
            cooler.merge_coolers(input_uris=merge, output_uri=merge_name, mergebuf=1000000)
        tmp_merge = cooler.Cooler(merge_name)
        cooler.balance_cooler(tmp_merge, store=True)
        resolution = res
        win = [10 * resolution]
        insulation_table = insulation(tmp_merge, win, verbose=False)
        boundary = insulation_table[insulation_table['is_boundary_' + str(win[0])] == True]
        
        boundary = boundary[boundary['chrom'].isin(chrom)]

        tad_index = list(boundary.index)

        result = []
        clr = self.clr_obj[self.resolution.index(res)]
        for i in tqdm.tqdm(range(len(clr)), desc='cell'):
            DS = self.__get_TAD_DS(clr[i], tad_index, epoch, nbins)
            result.append(DS)

        return result

    def __get_TAD_DS(self, clr, data, epoch, nbins):
        chrom = ['chr' + str(i) for i in range(20)]
        bins = clr.bins()[:]
        
        bins['index'] = list(bins.index)
        

        thre = bins.groupby('chrom').describe()['index'].loc[:, ['min', 'max']]#.astype(int)
        thre = thre[thre.index.isin(chrom)]
        thre['min'] = thre['min'].astype(int)
        thre['max'] = thre['max'].astype(int)
   
        pixels = clr.pixels()[:]
        pixels['chr1'] = list(bins.iloc[list(pixels['bin1_id']), 0])
        pixels['chr2'] = list(bins.iloc[list(pixels['bin2_id']), 0])
        pixels = pixels[pixels['chr1'] == pixels['chr2']]
        pixels = pixels[pixels['chr1'].isin(chrom)]
        pixels = pixels[pixels['chr2'].isin(chrom)]

        result = pd.DataFrame(columns=['tad_boundary', 'tad_left', 'tad_right', 'chrom'])
        result['tad_boundary'] = list(data)
        result['chrom'] = list(bins.loc[data, :]['chrom'])
        result = result.merge(thre, on='chrom', how='left')
        result['tad_left'] = [max(a, b - nbins) for a, b in zip(list(result['min']), data)]
        result['tad_right'] = [min(a, b + nbins) for a, b in zip(list(result['max']), data)]
        intra_pixels = pixels[pixels['bin1_id'] != pixels['bin2_id']]
        inter_pixels = pixels[pixels['bin1_id'] != pixels['bin2_id']]

        result['DS'] = result.apply(lambda x: self.__get_tad_ds(x, intra_pixels, inter_pixels, nbins), axis=1)

        thre['r_min'] = list(thre['min'] + nbins)
        thre['r_max'] = list(thre['max'] - nbins)

        rds = self.__get_tad_rds(thre, intra_pixels, inter_pixels, nbins, epoch)

        result = np.mean(result['DS']) - rds

        return result

    def __get_tad_rds(self, df, intra_pixels, inter_pixels, nbins, epoch):
        
        res = []
        for i in range(epoch):
            # print(nbins)
            idx = np.random.randint(df['r_min'], df['r_max'])
            idx = random.sample(list(idx), 1)[0]
            idx1 = range(idx - nbins, idx)
            idx2 = range(idx + 1, idx + nbins + 1)
            intra = intra_pixels[((intra_pixels['bin1_id'].isin(idx1)) & (intra_pixels['bin2_id'].isin(idx1))) | (
                        (intra_pixels['bin1_id'].isin(idx2)) & (intra_pixels['bin2_id'].isin(idx2)))]['count'].sum()
            inter = inter_pixels[(inter_pixels['bin1_id'].isin(idx1)) & (inter_pixels['bin2_id'].isin(idx2)) | (
                inter_pixels['bin1_id'].isin(idx2)) & (inter_pixels['bin2_id'].isin(idx1))]['count'].sum()
            # print(idx1,idx2)
            res.append((intra - inter) / (nbins * nbins * 2))

        return np.mean(res)

    def __get_tad_ds(self, x, intra_pixels, inter_pixels, nbins):
        x = list(x)

        idx1 = range(x[1], x[0])
        idx2 = range(x[0] + 1, x[2] + 1)

        intra = intra_pixels[((intra_pixels['bin1_id'].isin(idx1)) & (intra_pixels['bin2_id'].isin(idx1))) | (
                    (intra_pixels['bin1_id'].isin(idx2)) & (intra_pixels['bin2_id'].isin(idx2)))]['count'].sum()
        inter = inter_pixels[(inter_pixels['bin1_id'].isin(idx1)) & (inter_pixels['bin2_id'].isin(idx2)) | (
            inter_pixels['bin1_id'].isin(idx2)) & (inter_pixels['bin2_id'].isin(idx1))]['count'].sum()
        # print(intra, inter, 1)
        return (intra - inter) / (nbins * nbins * 2)

    def get_Centromere_nDs(self, res=None, epoch=100):
        if res is None:
            sys.exit('resoluiton is required')

        clr = self.clr_obj[self.resolution.index(res)]

        result = []
        for i, cr in enumerate(clr):
            log = self.__get_Centromere_nDs(cr, res, epoch=epoch)
            result.append([log[0], log[1], log[2]])
        result = pd.DataFrame(result)
        result.columns = ['DS', 'rDS', 'nDS']
        result.to_csv('./file/%s_centromere_DS_%s_%s.csv' % (self.describe, str(res), self.mcool_names[i]), sep='\t')

        return list(result['nDS'])

    def __get_Centromere_nDs(self, clr, res, epoch):
        centro_pos1 = 3000000  # 3mb
        centro_pos2 = 13000000  # 13mb

        # 获得binid的最大值和最小值，每个染色体不一定有这么长的bp区间
        bin_idx1 = int(centro_pos1 / res)
        bins_idx2 = int(centro_pos2 / res)

        bins = clr.bins()[:]
        pixels = clr.pixels()[:]

        chrom, n_bins = list(bins.groupby('chrom').count()['start'].index), list(bins.groupby('chrom').count()['start'])
        interval_to_cal = {}
        for c, n in zip(chrom, n_bins):
            if bin_idx1 > int(n):
                continue
            interval_to_cal[c] = [bin_idx1, min(n, bins_idx2)]

        ac = 0
        for i, k in enumerate(interval_to_cal.keys()):
            if i == 0:
                ac += n_bins[i]
            else:
                id1, id2 = interval_to_cal[k]
                id1 += ac
                id2 += ac
                ac += n_bins[i]
                interval_to_cal[k] = [id1, id2]
        interval_to_cal_ = []
        for k in interval_to_cal.keys():
            interval_to_cal_ += list(range(interval_to_cal[k][0], interval_to_cal[k][1] + 1))

        chr1 = interval_to_cal_
        chr2 = list(set(list(bins.index)).difference(set(chr1)))

        DS = self.__terri_detection_score(pixels, chr1, chr2)

        r_DS = []
        for i in range(epoch):
            random_interval = list(bins.sample(len(chr1)).index)
            random_interval2 = list(set(list(bins.index)).difference(set(random_interval)))
            r_DS.append(self.__terri_detection_score(pixels, random_interval, random_interval2))
        rDS = np.mean(r_DS)
        print(DS - rDS)

        return DS, rDS, DS - rDS

    def get_insulation_score(self, res=None, win=None):
        if res is None:
            sys.exit('resoluiton is required')

        clr = self.clr_obj[self.resolution.index(res)]
        result = []
        for i, cr in enumerate(clr):
            score = self.__get_insulation_score(cr, win)

        return result

    def __get_insulation_score(self, clr, win):
        win = [win]
        if not cooltools.lib.is_cooler_balanced(clr):
            cooler.balance_cooler(clr, store=True)
        # 获得 insulation score，以及是不是边界
        insulation_table = insulation(clr, win, verbose=True)
        s = insulation_table['log2_insulation_score_{}'.format(win[0])]
        s = np.nan_to_num(s)
        return s

    def pooled_cell(self, pool_size=None, num=None, sort=True):
        if self.__merge:
            sys.exit('merged can not pool the cell')
        if self.__pool:
            print('already pooled')
            return

        if pool_size is None or num is None:
            sys.exit('pool_size and num is required')

        # 一共多少个细胞
        element = len(self.clr_obj[0])
        #
        # if sort and num > element:
        #     print('sorted mode: num is bigger than element, set num = element')
        #     num = element
        new_clr_obj = []
        # 循环添加不同分辨率的空list
        for i in range(len(self.resolution)):
            new_clr_obj.append([])

        if not sort:
            for i in tqdm.tqdm(range(num), desc='pooling'):
                pool = random.choices(list(range(element)), k=pool_size)
                mcool_pool = [self.cache_file[idx] for idx in pool]
                for j in range(len(self.resolution)):
                    clr_list = [val + '::/resolutions/' + str(self.resolution[j]) for val in mcool_pool]

                    new_clr_obj[j].append(self.__merge_pool(clr_list, str(i)))
        else:
            # 从大到小排序total contact
            idx = [i[0] for i in sorted(enumerate(self.get_total_contact()), key=lambda x: x[1], reverse=True)]
            # 如果num大于element，那么就随机插入
            if num > element:
                # 假设这是原始列表
                original_list = self.get_total_contact()
                selected_elements = original_list.copy()
                selected_elements.extend(random.sample(original_list, 5))
                sorted_elements = sorted(selected_elements, reverse=True)
                idx = [original_list.index(x) for x in sorted_elements]

            for i in tqdm.tqdm(range(num), desc='pooling::sorted'):
                pool = idx[:i + 1]

                mcool_pool = [self.cache_file[idx] for idx in pool]
                for j in range(len(self.resolution)):
                    clr_list = [val + '::/resolutions/' + str(self.resolution[j]) for val in mcool_pool]
                    new_clr_obj[j].append(self.__merge_pool(clr_list, str(i)))

        self.clr_obj = new_clr_obj
        self.mcool_names = ['pooled_cell_%s' % str(i) for i in range(num)]

    def __merge_pool(self, clr_list, cell_ids):

        cooler.merge_coolers(input_uris=clr_list,
                             output_uri='./cache/%s_pool_all_res-%s.mcool' % (self.describe, cell_ids),
                             mergebuf=1000000)

        res = cooler.Cooler('./cache/%s_pool_all_res-%s.mcool' % (self.describe, cell_ids))

        return res

    def filter_by_total_contact(self, max_=np.inf, min_=0):
        if max_ is None and min_ is None:
            sys.exit('max and min is required')

        tmp_mcool_files = []
        tmp_mcool_names = []
        tmp_cache_file = []
        tmp_clr_obj = []
        for i in range(len(self.resolution)):
            tmp_clr_obj.append([])
        total_contact_list = self.get_total_contact()
        for a, b, c, in zip(self.mcool_files, self.mcool_names, self.cache_file):
            # print(c)
            total_contact = total_contact_list[self.mcool_files.index(a)]
            if total_contact > max_ or total_contact < min_:
                continue
            tmp_mcool_files.append(a)
            tmp_mcool_names.append(b)
            tmp_cache_file.append(c)
            for i in range(len(self.resolution)):
                tmp_clr_obj[i].append(self.clr_obj[i][self.mcool_files.index(a)])

        self.mcool_files = tmp_mcool_files
        self.mcool_names = tmp_mcool_names
        self.cache_file = tmp_cache_file
        self.clr_obj = tmp_clr_obj

    def get_cis_trans_ratio(self):


        result = []


        for clr in self.clr_obj[0]:
            cis, trans = self.__get_cis_trans_ratio(clr)
            result.append([cis, trans])

        result = pd.DataFrame(result)
        result.columns = ['cis', 'trans']
        result.index = self.mcool_names
        return result

    def __get_cis_trans_ratio(self, clr):
        bins = clr.bins()[:]
        pixels = clr.pixels()[:]
        pixels['chr1'] = list(bins.iloc[list(pixels['bin1_id']), 0])
        pixels['chr2'] = list(bins.iloc[list(pixels['bin2_id']), 0])

        cis = sum(pixels[pixels['chr1'] == pixels['chr2']]['count'])
        trans = sum(pixels[pixels['chr1'] != pixels['chr2']]['count'])

        return cis, trans

    def sorted_cell(self, reverse=True, num=None):
        idx = [i[0] for i in sorted(enumerate(self.get_total_contact()), key=lambda x: x[1], reverse=reverse)]

        self.cache_file = [self.cache_file[i] for i in idx]
        self.mcool_files = [self.mcool_files[i] for i in idx]
        self.mcool_names = [self.mcool_names[i] for i in idx]
        for i in range(len(self.resolution)):
            self.clr_obj[i] = [self.clr_obj[i][j] for j in idx]

        if num is not None and num < len(self.cache_file):
            self.cache_file = self.cache_file[:num]
            self.mcool_files = self.mcool_files[:num]
            self.mcool_names = self.mcool_names[:num]
            for i in range(len(self.resolution)):
                self.clr_obj[i] = self.clr_obj[i][:num]


    # plot coverage rate
    def plot_coverage(self, res=None, save=False):

        if res is None:
            df = []
            for i,res_clr in enumerate(self.clr_obj):
                for clr in res_clr:
                    cis, total = cooltools.coverage(clr, ignore_diags=1)
                    coverage = sum(cis > 0) / cis.shape[0]
                    df.append([coverage, self.resolution[i]])
            df = pd.DataFrame(df)
            df.columns = ['coverage', 'resolution']


        else:
            if res not in self.resolution:
                sys.exit('res is not in resolution')
            idx = self.resolution.index(res)

            res_clr = self.clr_obj[idx]
            df = []
            for clr in res_clr:
                cis, total = cooltools.coverage(clr, ignore_diags=1)
                coverage = sum(cis > 0) / cis.shape[0]
                df.append([coverage, res])
            df = pd.DataFrame(df)
            df.columns = ['coverage', 'resolution']

        fig = plt.figure(figsize=(10, 6),dpi=100)
        sns.violinplot(x='resolution', y='coverage', data=df)

        if save:
            if isinstance(save, bool):
                os.makedirs('./figures', exist_ok=True)
                plt.savefig('./figures/violinplot_coverage.png')
            elif isinstance(save, str):
                plt.savefig('./figures/%s_violinplot_coverage.png' % save)
            else:
                sys.exit('save must be bool or str')


        plt.show()



        return df

    # plot short-mid-long
    def plot_range(self,  mode='long', save=False):
        '''


        Parameters
        ----------
        res
        mode
            short : require resolution <= 10000, which is < 10kb
            mid : require resolution > 10000 and <= 100000, which is > 20kb and <= 2Mb
            long : require resolution > 100000, which is > 2Mb

        save

        Returns
        -------

        '''
        if 10000 not in self.resolution:
            sys.exit('resolution must contain 10000')
        df = []
        for i,clr in enumerate(self.clr_obj[0]):
            bins = clr.bins()[:]
            pixels = clr.pixels()[:]
            pixels['chr1'] = list(bins.iloc[list(pixels['bin1_id']), 0])
            pixels['chr2'] = list(bins.iloc[list(pixels['bin2_id']), 0])
            pixels = pixels[pixels['chr1'] == pixels['chr2']]
            total = clr.info['sum']


            pixels['bin1_id'] = pixels['bin1_id'].astype(int)
            pixels['bin2_id'] = pixels['bin2_id'].astype(int)

            short = sum(pixels[pixels['bin2_id'] - pixels['bin1_id'] <= 2]['count'])

            mid = sum(
                pixels[(pixels['bin2_id'] - pixels['bin1_id'] > 2) & (pixels['bin2_id'] - pixels['bin1_id'] <= 200)][
                    'count'])

            long = sum(pixels[pixels['bin2_id'] - pixels['bin1_id'] > 200]['count'])
            df.append([short / total, self.mcool_names[i], 'short'])
            df.append([mid / total, self.mcool_names[i], 'mid'])
            df.append([long / total, self.mcool_names[i], 'long'])

        df = pd.DataFrame(df)
        df.columns = ['ratio', 'name', 'range']


        fig = plt.figure(figsize=(10, 6),dpi=100)
        sns.violinplot(x='range', y='ratio', data=df)

        if save:
            if isinstance(save, bool):
                os.makedirs('./figures', exist_ok=True)
                plt.savefig('./figures/violinplot_range.png')
            elif isinstance(save, str):
                plt.savefig('./figures/%s_violinplot_range.png' % save)
            else:
                sys.exit('save must be bool or str')


        return df


    def save(self, path, worker=2):

        os.makedirs(path, exist_ok=True)
        os.makedirs('./cache', exist_ok=True)

        task = []
        for clr,name in zip(self.clr_obj[0], self.mcool_names):
            task.append([clr,name,path])

        Parallel(n_jobs=worker)(delayed(self.save__)(clr,name,path,worker) for clr,name,path in task)


    def save__(self, clr, name, path,worker):
        cooler.create_cooler(
            './cache/%s.cool' % name,
            bins=clr.bins()[:],
            pixels=clr.pixels()[:]
        )

        cooler.zoomify_cooler(base_uris='./cache/%s.cool' % name, outfile=path + '/' + name + '.mcool',
                              resolutions=self.resolution,
                              chunksize=10000, nproc=worker)
        os.remove('./cache/%s.cool' % name)


    # QC
    def filter_by_coverage(self, min_=0, max_=1):
        if min_ > max_:
            sys.exit('min must be smaller than max')

        tmp_mcool_files = []
        tmp_mcool_names = []
        tmp_cache_file = []
        tmp_clr_obj = []
        for i in range(len(self.resolution)):
            tmp_clr_obj.append([])
        for a, b, c, in zip(self.mcool_files, self.mcool_names, self.cache_file):
            cis, total = cooltools.coverage(a + '::/resolutions/' + str(self.resolution[0]), ignore_diags=1)
            coverage = sum(cis > 0) / cis.shape[0]
            if coverage < min_ or coverage > max_:
                continue
            tmp_mcool_files.append(a)
            tmp_mcool_names.append(b)
            tmp_cache_file.append(c)
            for i in range(len(self.resolution)):
                tmp_clr_obj[i].append(self.clr_obj[i][self.mcool_files.index(a)])

        self.mcool_files = tmp_mcool_files
        self.mcool_names = tmp_mcool_names
        self.cache_file = tmp_cache_file
        self.clr_obj = tmp_clr_obj

    def filter_by_range(self, s_min_=0, s_max_=1,
                                m_min_=0, m_max_=1,
                                l_min_=0, l_max_=1):
        if 10000 not in self.resolution:
            sys.exit('resolution must contain 10000')
        if s_min_ > s_max_ or m_min_ > m_max_ or l_min_ > l_max_:
            sys.exit('min must be smaller than max')

        tmp_mcool_files = []
        tmp_mcool_names = []
        tmp_cache_file = []
        tmp_clr_obj = []
        for i in range(len(self.resolution)):
            tmp_clr_obj.append([])
        for a, b, c, in zip(self.mcool_files, self.mcool_names, self.cache_file):
            clr = cooler.Cooler(c)
            bins = clr.bins()[:]
            pixels = clr.pixels()[:]
            pixels['chr1'] = list(bins.iloc[list(pixels['bin1_id']), 0])
            pixels['chr2'] = list(bins.iloc[list(pixels['bin2_id']), 0])
            pixels = pixels[pixels['chr1'] == pixels['chr2']]
            total = clr.info['sum']

            pixels['bin1_id'] = pixels['bin1_id'].astype(int)
            pixels['bin2_id'] = pixels['bin2_id'].astype(int)

            short = sum(pixels[pixels['bin2_id'] - pixels['bin1_id'] <= 2]['count'])

            mid = sum(
                pixels[(pixels['bin2_id'] - pixels['bin1_id'] > 2) & (pixels['bin2_id'] - pixels['bin1_id'] <= 200)][
                    'count'])

            long = sum(pixels[pixels['bin2_id'] - pixels['bin1_id'] > 200]['count'])
            l_ratio = long / total
            m_ratio = mid / total
            s_ratio = short / total
            if l_ratio < l_min_ or l_ratio > l_max_:
                continue

            if m_ratio < m_min_ or m_ratio > m_max_:
                continue

            if s_ratio < s_min_ or s_ratio > s_max_:
                continue

            tmp_mcool_files.append(a)
            tmp_mcool_names.append(b)
            tmp_cache_file.append(c)
            for i in range(len(self.resolution)):
                tmp_clr_obj[i].append(self.clr_obj[i][self.mcool_files.index(a)])

        self.mcool_files = tmp_mcool_files
        self.mcool_names = tmp_mcool_names
        self.cache_file = tmp_cache_file


    def filter_by_cistrans(self, min_=0, max_=np.inf):
        if min_ > max_:
            sys.exit('min must be smaller than max')

        tmp_mcool_files = []
        tmp_mcool_names = []
        tmp_cache_file = []
        tmp_clr_obj = []
        for i in range(len(self.resolution)):
            tmp_clr_obj.append([])
        for a, b, c, in zip(self.mcool_files, self.mcool_names, self.cache_file):
            clr = cooler.Cooler(a + '::/resolutions/' + str(self.resolution[0]))
            cis, trans = self.__get_cis_trans_ratio(clr)
            ratio = cis/trans
            if ratio < min_ or ratio > max_:
                continue

            tmp_mcool_files.append(a)
            tmp_mcool_names.append(b)
            tmp_cache_file.append(c)
            for i in range(len(self.resolution)):
                tmp_clr_obj[i].append(self.clr_obj[i][self.mcool_files.index(a)])

        self.mcool_files = tmp_mcool_files
        self.mcool_names = tmp_mcool_names
        self.cache_file = tmp_cache_file


    def get_giniQC(self, n_jobs):
        print('balancing')
        
        for clr in self.clr_obj[-1]:
            if not cooltools.lib.is_cooler_balanced(clr):
                cooler.balance_cooler(clr, store=True)
        result = []
        for clr in self.clr_obj[-1]:
            tmp = calculate_giniqc(clr)
            result.append(tmp)

        # result = Parallel(n_jobs=n_jobs)(calculate_giniqc(clr) for clr in self.clr_obj[-1])
        print('balanced')
        result = pd.DataFrame(result)
        print(result)
        result.columns = ['gini', 'adjusted_gini']

        result.index = self.mcool_names

        return result





   

































