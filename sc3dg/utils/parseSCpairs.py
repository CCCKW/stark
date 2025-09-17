#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import gzip

def open_file(filename, mode='rt'):
    """根据文件扩展名自动选择打开方式"""
    if filename.endswith('.gz'):
        return gzip.open(filename, mode=mode, encoding='utf-8')
    else:
        return open(filename, mode=mode, encoding='utf-8')

def main(input_file, output_dir="./scPairDir"):
    """Split .pairs or .pairs.gz file by CellBarcode (extracted from 2nd and 3rd underscore-separated fields)."""
    
    # 初始化变量
    Head_info = []
    scPairDict = {}

    # 读取输入文件
    try:
        with open_file(input_file) as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line:
                    continue  # 跳过空行
                if line.startswith("#"):
                    Head_info.append(line + "\n")
                else:
                    parts = line.split("\t")[0].split("_")
                    if len(parts) < 3:
                        print(f"Warning: Line {line_num} skipped (not enough parts): {line}")
                        continue
                    CellBarcodes = "_".join(parts[1:])
                    #print(CellBarcodes)
                    if CellBarcodes not in scPairDict:
                        scPairDict[CellBarcodes] = []
                    scPairDict[CellBarcodes].append(line + "\n")
    except FileNotFoundError:
        print(f"Error: File '{input_file}' not found.")
        exit(1)
    except Exception as e:
        print(f"Error reading file '{input_file}': {e}")
        exit(1)

    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)

    # 写入每个 barcode 的文件
    for barcode, pairinfo in scPairDict.items():
        output_path = os.path.join(output_dir, f"{barcode}.pairs")
        try:
            with open(output_path, 'w', encoding='utf-8') as fp:
                fp.writelines(Head_info)
                fp.writelines(pairinfo)
        except Exception as e:
            print(f"Error writing file {output_path}: {e}")

    print(f"Success! Split {len(scPairDict)} barcode groups into {output_dir}/")

if __name__ == "__main__":
    # 定义参数变量
    input_file = "input.pairs.gz"  # 替换为实际的输入文件路径
    output_dir = "./scPairDir"     # 替换为实际的输出目录路径
    
    main(input_file, output_dir)