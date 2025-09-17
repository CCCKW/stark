

import click
from joblib import Parallel, delayed
import cooler
import os
import sc3dg.cooltools.cooltools as cooltools
import pandas as pd
from sc3dg.analysis.api import MMCooler
from sc3dg.cooltools.cooltools import insulation
chrom = ['chr' + str(i) for i in range(20)]
@click.command('accum', short_help='accum analysis')
@click.option('--mcool',  help='mcool path',required=True, default=None)
@click.option('--output',  help='output path',required=True, default=None)
@click.option('--resolution', help='resolution',required=True, default=40000, type=int)
@click.option('--top_num', help='top number',required=True, default=5, type=int)
def accum(mcool, output, resolution, top_num, describe='accum'):
    pass
    ob1 = MMCooler(directory=mcool,
                   resolution=[resolution],
                   describe = describe,
                   merge=False
                   )
    ob1.sorted_cell(num=top_num)

    iou0 = []  # for different windowns
    iou1 = []
    iou2 = []
    iou3 = []
    iou = []
    win = [3 * resolution, 5 * resolution, 10 * resolution, 25 * resolution]
    print('merging')
    cache_file = []
    for clr in ob1.mcool_files:
        cache_file.append(clr + '::resolutions/%s' % resolution)

    if not os.path.exists('./cache/%s_assemb.mcool' % (describe)):
        print('merging')
        cooler.merge_coolers(input_uris=cache_file, output_uri='./cache/%s_assemb.mcool' % (describe),
                             mergebuf=10000000000)

    clr = cooler.Cooler('./cache/%s_assemb.mcool' % (describe))
    if not cooltools.lib.is_cooler_balanced(clr):
        cooler.balance_cooler(clr, mad_max=0, min_nnz=0, ignore_diags=1
                              , store=True, store_name='weight')

    insulation_table = insulation(clr, win, verbose=True, nproc=20)
    boundary = insulation_table[insulation_table['is_boundary_' + str(win[0])] == True]
    boundary2 = boundary[boundary['chrom'].isin(chrom)]

    print(boundary2)

    num = len(ob1.clr_obj[0])
    tmp0 = []
    tmp1 = []
    tmp2 = []
    tmp3 = []
    tmp = []
    count0 = []
    count1 = []
    count2 = []
    count3 = []



    ob1.pooled_cell(pool_size=1000, num=min(num, top_num))

    task = []
    for i, clr in enumerate(ob1.clr_obj[-1]):
        print(i)
        task.append((clr, i))

    def balance(clr, c):
        print(clr, c)
        if not cooltools.lib.is_cooler_balanced(clr):
            cooler.balance_cooler(clr, mad_max=0, min_nnz=0, ignore_diags=1
                                  , store=True, store_name='weight')

    Parallel(n_jobs=10)(delayed(balance)(clr, c) for clr, c in task)

    for i, clr in enumerate(ob1.clr_obj[-1]):
        print(i)
        insulation_table = insulation(clr, win, verbose=False, )
        boundary = insulation_table[insulation_table['is_boundary_' + str(win[0])] == True]
        boundary = boundary[boundary['chrom'].isin(chrom)]
        b0 = boundary[boundary['is_boundary_' + str(win[0])] == True]
        b1 = boundary[boundary['is_boundary_' + str(win[1])] == True]
        b2 = boundary[boundary['is_boundary_' + str(win[2])] == True]
        b3 = boundary[boundary['is_boundary_' + str(win[3])] == True]

        b00 = boundary2[boundary2['is_boundary_' + str(win[0])] == True]
        b11 = boundary2[boundary2['is_boundary_' + str(win[1])] == True]
        b22 = boundary2[boundary2['is_boundary_' + str(win[2])] == True]
        b33 = boundary2[boundary2['is_boundary_' + str(win[3])] == True]

        iou_0 = len(set(b0.index).intersection(set(b00.index))) / len(set(b0.index).union(set(b00.index)))
        iou_1 = len(set(b1.index).intersection(set(b11.index))) / len(set(b1.index).union(set(b11.index)))
        iou_2 = len(set(b2.index).intersection(set(b22.index))) / len(set(b2.index).union(set(b22.index)))
        iou_3 = len(set(b3.index).intersection(set(b33.index))) / len(set(b3.index).union(set(b33.index)))
        iou_ = len(set(boundary.index).intersection(set(boundary2.index))) / len(
            set(boundary.index).union(set(boundary2.index)))

        tmp0.append(iou_0)
        tmp1.append(iou_1)
        tmp2.append(iou_2)
        tmp3.append(iou_3)
        tmp.append(iou_)

        count0.append(b0.shape[0])
        count1.append(b1.shape[0])
        count2.append(b2.shape[0])
        count3.append(b3.shape[0])

    iou0.append(tmp0)
    iou1.append(tmp1)
    iou2.append(tmp2)
    iou3.append(tmp3)
    iou.append(tmp)

    res = pd.DataFrame(iou0)
    res.to_csv('%s/bulkvspooled_0_%s.csv' % (output, describe), index=False)
    res = pd.DataFrame(iou1)
    res.to_csv('%s/bulkvspooled_1_%s.csv' % (output, describe), index=False)
    res = pd.DataFrame(iou2)
    res.to_csv('%s/bulkvspooled_2_%s.csv' % (output, describe), index=False)
    res = pd.DataFrame(iou3)
    res.to_csv('%s/bulkvspooled_3_%s.csv' % (output, describe), index=False)


    res = pd.DataFrame(count0)
    res.to_csv('%s/bulkvspooled_0_%s_count.csv' % (output, describe), index=False)
    res = pd.DataFrame(count1)
    res.to_csv('%s/bulkvspooled_1_%s_count.csv' % (output, describe), index=False)
    res = pd.DataFrame(count2)
    res.to_csv('%s/bulkvspooled_2_%s_count.csv' % (output, describe), index=False)
    res = pd.DataFrame(count3)
    res.to_csv('%s/bulkvspooled_3_%s_count.csv' % (output, describe), index=False)
    res = pd.DataFrame(iou)
    res.to_csv('%s/bulkvspooled_%s.csv' % (output, describe), index=False)




