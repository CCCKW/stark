#!/usr/bin/env python3

import sys
import gzip

def split_pairs_by_barcode(pairfile, output_path):
    """
    将 .pairs 文件（支持 .gz 压缩）按 barcode 拆分为多个文件。
    输出记录数最多的前 2000 个 barcode，每个文件包含 header。

    Args:
        pairfile (str): 输入的 .pairs 文件路径（可为 .gz）
        output_path (str): 输出目录（需以 / 结尾）
    """
    PrirDict = {}          # barcode -> lines
    Head_info = []         # header 行（以 # 开头）

    print(f"🔍 正在读取文件: {pairfile}")

    # 判断是否为 gzip 压缩文件
    open_func = gzip.open if pairfile.endswith('.gz') else open
    mode = 'rt' if pairfile.endswith('.gz') else 'r'  # gzip 需要文本模式
    count = 0
    try:
        with open_func(pairfile, mode) as f:
            for line in f:
                count+=1
                if count %10000000==0:
                    print(f"has processed {count} reads")
                line = line.strip('\n')
                if line.startswith('#'):
                    Head_info.append(line + '\n')
                    continue

                if not line:
                    continue

                # 解析第一列：barcode:XXXXXX
                try:
                    barcode_field = line.split("\t")[0]
                    barcode = barcode_field.split(":", 1)[1]  # 取冒号后部分
                except (IndexError, ValueError):
                    print(f"⚠️ 跳过格式错误的行: {line}", file=sys.stderr)
                    continue

                if barcode not in PrirDict:
                    PrirDict[barcode] = []
                PrirDict[barcode].append(line + '\n')

        print(f"✅ 共读取到 {len(PrirDict)} 个唯一 barcode")

        # 按记录数量排序，取前 2000 名
        sorted_barcodes = sorted(PrirDict.items(), key=lambda x: len(x[1]), reverse=True)
        top_Cell = sorted_barcodes[:10000]

        print(f"🏆 选取记录数最多的前 2000 个 barcode")

        # 写出每个选中的 barcode 文件
        for barcode, lines in top_Cell:
            output_file = f"{output_path}{barcode}.pairs"
            try:
                with open(output_file, 'w') as fp:
                    fp.writelines(Head_info)
                    fp.writelines(lines)
                print(f"📄 已生成文件: {output_file} ({len(lines)} 条记录)")
            except Exception as e:
                print(f"❌ 写入文件失败 {output_file}: {e}", file=sys.stderr)

    except FileNotFoundError:
        print(f"❌ 错误：找不到文件 '{pairfile}'", file=sys.stderr)
        sys.exit(1)
    except OSError as e:
        print(f"❌ 读取文件时发生 I/O 错误: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"❌ 处理过程中发生未知错误: {e}", file=sys.stderr)
        sys.exit(1)


# ============ 主程序 ============
if __name__ == '__main__':
    # 检查参数
    if len(sys.argv) != 3:
        print("📌 用法: python split_pairs.py <input.pairs[.gz]> <output_dir/>", file=sys.stderr)
        print("示例:", file=sys.stderr)
        print("    python split_pairs.py 1.pairs ./split/", file=sys.stderr)
        print("    python split_pairs.py 1.pairs.gz ./split/", file=sys.stderr)
        sys.exit(1)

    pairfile = sys.argv[1]
    output_path = sys.argv[2]

    # 确保 output_path 以 / 结尾
    if not output_path.endswith('/'):
        output_path += '/'

    split_pairs_by_barcode(pairfile, output_path)
