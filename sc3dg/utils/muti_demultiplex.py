import gzip, os, sys, datetime, itertools, pickle, re
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed



# ================== 全局常量 ==================
hash_base = 5
hash_power = hash_base ** np.arange(10)

# ================== 工具函数 ==================
def HASH(b):
    return (b % hash_base) @ hash_power[:len(b)]

def open_file(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

def print_datetime():
    return datetime.datetime.now().strftime('[%Y/%m/%d %H:%M:%S]\t')

def matchBarcode(seq, bc, mode):
    if isinstance(seq, bytes):
        seq = np.frombuffer(seq, dtype=np.uint8)
    if bc.ndim == 1:
        return bc[HASH(seq)]
    else:
        d = ((seq[None] != bc) & (bc != ord('N'))).sum(1)
        idx = np.argmin(d)
        if (mode == 'unique' or d[idx] <= mode) and (d == d[idx]).sum() == 1:
            return idx
        else:
            return -1

def preprocessBarcode(barcode, mode):
    N, L = barcode.shape
    d = np.full(hash_base**L, -1)
    if mode == 0:
        for i, b in enumerate(barcode):
            d[HASH(b)] = i
    elif mode == 1:
        for i, b in enumerate(barcode):
            bb = np.copy(b)
            for idx in range(L):
                for c in [65, 67, 71, 78, 84]:  # A, C, G, N, T
                    bb[idx] = c
                    d[HASH(bb)] = i
                bb[idx] = b[idx]
    elif mode == 'unique' or mode == 2:
        b = np.zeros(L, dtype=np.int64)
        def f(i):
            if i == L:
                d[HASH(b)] = matchBarcode(b, barcode, mode)
            else:
                for c in [65, 67, 71, 78, 84]:
                    b[i] = c
                    f(i + 1)
        f(0)
    else:
        raise NotImplementedError
    return d

def loadBarcode(filename, mode, preprocessed_barcode_file_prefix=None):
    with open(filename, 'r') as f:
        ret = [line.rstrip('\n').split('\t') for line in f]
    N = len(ret)
    wellID, barcode1, barcode2 = zip(*ret)
    barcode1 = np.array([np.frombuffer(bytes(_, 'utf-8'), dtype=np.uint8) for _ in barcode1])
    barcode2 = np.array([np.frombuffer(bytes(_, 'utf-8'), dtype=np.uint8) for _ in barcode2])
    barcodes = []
    for i, bc in [(1, barcode1), (2, barcode2)]:
        while True:
            if preprocessed_barcode_file_prefix is not None:
                filename = preprocessed_barcode_file_prefix + f'{i}.pkl'
            else:
                filename = None
            if filename is not None and os.path.exists(filename):
                sys.stderr.write(f'{print_datetime()}Loading preprocessed barcodes from {filename}\n')
                with open(filename, 'rb') as f:
                    bc = pickle.load(f)
                break
            if N > 2:
                bc = preprocessBarcode(bc, mode)
                if filename is not None:
                    with open(filename, 'wb') as f:
                        pickle.dump(bc, f)
                    sys.stderr.write(f'{print_datetime()}Saving preprocessed barcodes to {filename}\n')
                break
            else:
                break
        barcodes.append(bc)
        if bc.ndim == 1:
            sys.stderr.write(f'{print_datetime()}R{i}: {(bc != -1).sum()} out of {len(bc)} seq are valid\n')
    barcode1, barcode2 = barcodes
    return N, wellID, barcode1, barcode2

def is_slice_empty(s):
    if s is None:
        return True
    start, stop, step = s.indices(1000000)  # 用一个大数模拟长度
    return start >= stop

# ================== 多进程处理单元 ==================
def process_chunk(args):
    """
    处理一个 chunk 的 read1 和 read2 数据
    返回: (results, num_matched)
    """
    chunk, wellID_list, barcode1_list, barcode2_list, pos1_list, pos2_list, mode_list, barcode2keep = args

    K = len(mode_list)
    num_matched = np.zeros([K+1, 3], dtype=np.uint64)  # R1, R2, both
    results = []  # 存储要写入的数据 (out1, out2, out1_u, out2_u)

    for lines_read1, lines_read2 in chunk:
        flag_list = []
        idx_list = []

        for wellID, bc1, bc2, pos1, pos2, mode in zip(wellID_list, barcode1_list, barcode2_list,
                                                      pos1_list, pos2_list, mode_list):
            seq1 = np.frombuffer(lines_read1[1], dtype=np.uint8)[pos1] if len(pos1) > 0 else None
            if not is_slice_empty(pos2):
                seq2 = np.frombuffer(lines_read2[1], dtype=np.uint8)[pos2]
            else:
                seq2 = None

            idx1 = matchBarcode(seq1, bc1, mode) if seq1 is not None else None
            idx2 = matchBarcode(seq2, bc2, mode) if seq2 is not None else None

            matched_r1 = (idx1 != -1) if idx1 is not None else True
            matched_r2 = (idx2 != -1) if idx2 is not None else True
            matched_both = (matched_r1 and matched_r2 and
                            (idx1 == idx2 or idx1 is None or idx2 is None))

            flag_list.append((matched_r1, matched_r2, matched_both))
            if matched_both and idx_list is not None:
                idx = idx1 if idx1 is not None else (idx2 if idx2 is not None else -1)
            else:
                idx = -1
            idx_list.append(idx)

        flag_list = np.array(flag_list, dtype=bool)
        flag_list_all = flag_list.all(axis=0)  # AND across all barcodes

        if flag_list_all[-1]:  # 所有 barcode 都匹配
            wellID = '_'.join(wellID_list[i][idx_list[i]] for i in range(K))
            if barcode2keep is not None and wellID not in barcode2keep:
                flag_list_all[-1] = False

        if flag_list_all[-1]:
            t = ('_' + wellID + '\n').encode('utf-8')
            out1 = b''.join([lines_read1[0].rstrip(b'\n').split(b' ')[0] + t] + lines_read1[1:])
            out2 = b''.join([lines_read2[0].rstrip(b'\n').split(b' ')[0] + t] + lines_read2[1:])
            out1_u = out2_u = None
        else:
            out1 = out2 = b'\n' * 4
            out1_u = b''.join(lines_read1)
            out2_u = b''.join(lines_read2)

        results.append((out1, out2, out1_u, out2_u))
        num_matched[:-1] += flag_list
        num_matched[-1] += flag_list_all

    return results, num_matched

# ================== 主函数 ==================
def main(infile1, infile2, outfile1, outfile2, barcode, pos1, pos2, mode, unknown1=None, unknown2=None, barcode2keep=None, chunk_size=100000, nproc=4):
    sys.stderr.write(f'{print_datetime()}Begin\n')
    sys.stderr.write(f'{print_datetime()}Parameters: infile1={infile1}, infile2={infile2}, outfile1={outfile1}, outfile2={outfile2}\n')

    # 加载 barcode
    wellID_list = []
    barcode1_list = []
    barcode2_list = []
    for barcode_file, mode_val in zip(barcode, mode):
        n, wellID, bc1, bc2 = loadBarcode(barcode_file, mode_val)
        wellID_list.append(wellID)
        barcode1_list.append(bc1)
        barcode2_list.append(bc2)
        sys.stderr.write(f'{print_datetime()}Loaded {n} barcodes from {barcode_file}\n')

    # 加载 barcode2keep
    if barcode2keep:
        with open(barcode2keep, 'r') as f:
            barcode2keep = {line.strip().replace('\t', '_') for line in f}
        sys.stderr.write(f'{print_datetime()}Keep {len(barcode2keep)} barcodes\n')
    else:
        barcode2keep = None
        sys.stderr.write(f'{print_datetime()}Keep all barcodes\n')

    # 打开文件
    in1 = open_file(infile1, 'rb')
    in2 = open_file(infile2, 'rb')
    out1 = open_file(outfile1, 'wb')
    out2 = open_file(outfile2, 'wb')

    if unknown1 and unknown2:
        out_u1 = open_file(unknown1, 'wb')
        out_u2 = open_file(unknown2, 'wb')
    else:
        out_u1 = out_u2 = None

    # 初始化统计
    total_matched = np.zeros([len(mode) + 1, 3], dtype=np.uint64)
    total_count = 0

    # 多进程池
    with ProcessPoolExecutor(max_workers=nproc) as executor:
        futures = []

        while True:
            chunk = []
            try:
                for _ in range(chunk_size):
                    lines1 = [in1.readline() for _ in range(4)]
                    lines2 = [in2.readline() for _ in range(4)]
                    if not lines1[0] or not lines2[0]:  # EOF
                        break
                    chunk.append((lines1, lines2))
                if not chunk:
                    break

                # 提交任务
                future = executor.submit(
                    process_chunk,
                    (chunk, wellID_list, barcode1_list, barcode2_list,
                     pos1, pos2, mode, barcode2keep)
                )
                futures.append(future)

                # 当任务数达到进程数时，开始回收结果
                if len(futures) >= nproc:
                    for future in as_completed(futures):
                        results, stats = future.result()
                        futures.remove(future)

                        # 写入结果
                        for out1_data, out2_data, out1_u, out2_u in results:
                            out1.write(out1_data)
                            out2.write(out2_data)
                            if out_u1 and out1_u:
                                out_u1.write(out1_u)
                                out_u2.write(out2_u)

                        # 更新统计
                        total_matched += stats
                        total_count += len(results)

                        if total_count % 10_000_000 == 0:
                            sys.stderr.write(f'{print_datetime()}# = {total_count}\n')
                            rate = total_matched.T / total_count
                            output = np.array2string(rate, separator='\t', formatter={'all': '{:.4f}'.format})
                            # 手动替换掉不需要的字符
                            output = output.replace('[', '').replace(']', '').replace(' ', '').replace('\n', '')
                            sys.stderr.write(output + '\n')
                            sys.stderr.flush()

            except Exception as e:
                sys.stderr.write(f'{print_datetime()}Error: {e}\n')
                break

        # 收集剩余任务
        for future in as_completed(futures):
            results, stats = future.result()
            for out1_data, out2_data, out1_u, out2_u in results:
                out1.write(out1_data)
                out2.write(out2_data)
                if out_u1 and out1_u:
                    out_u1.write(out1_u)
                    out_u2.write(out2_u)
            total_matched += stats
            total_count += len(results)

    # 关闭文件
    for f in [in1, in2, out1, out2]:
        if hasattr(f, 'close'): f.close()
    if out_u1: out_u1.close()
    if out_u2: out_u2.close()

    # 正确写法：使用 replace 代替跨行正则
    sys.stderr.write(f'{print_datetime()}Total # = {total_count}\n')
    sys.stderr.write(f'{print_datetime()}Final match stats:\n')

    # 先转成字符串，再清理
    output = np.array2string(total_matched.T, separator='\t')
    output = output.replace('[', '').replace(']', '').replace(' ', '')
    sys.stderr.write(output + '\n')

    sys.stderr.write(f'{print_datetime()}Done\n')

if __name__ == '__main__':
    # 定义参数变量 - 对应您的原始命令行参数
    # 您需要将这些路径替换为实际的文件路径
    adaptor_file_bc2 = "path/to/bc2.txt"  # 替换为实际路径
    adaptor_file_l2 = "path/to/l2.txt"    # 替换为实际路径
    adaptor_file_bc1 = "path/to/bc1.txt"  # 替换为实际路径
    fq_r1 = "sample_R1.fastq.gz"          # 替换为实际文件名
    fq_r2 = "sample_R2.fastq.gz"          # 替换为实际文件名
    faq = "path/to/fastq/dir"             # 替换为实际目录路径
    
    infile1 = f'{faq}/{fq_r1}'
    infile2 = f'{faq}/{fq_r2}'
    outfile1 = f'{faq}/demulti_R1.fastq.gz'
    outfile2 = f'{faq}/demulti_R2.fastq.gz'
    barcode = [adaptor_file_bc2, adaptor_file_l2, adaptor_file_bc1]
    pos1 = [[], [], []]  # 对应 "[]" 参数
    pos2 = [slice(0, 8), slice(8, 23), slice(23, 31)]
    mode = [1, 5, 1]
    unknown1 = None
    unknown2 = None
    barcode2keep = None  # 对应 --barcode2keep=""
    chunk_size = 100000
    nproc = 8
    
    main(infile1, infile2, outfile1, outfile2, barcode, pos1, pos2, mode, 
         unknown1, unknown2, barcode2keep, chunk_size, nproc)