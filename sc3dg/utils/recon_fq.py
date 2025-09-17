#!/usr/bin/env python3
"""
sam_to_fastq.py

将 SAM 文件转换为两个 FASTQ 文件（R1 和 R2），按行顺序交替写入。
注意：此方法假设配对 reads 在文件中成对且顺序连续（.1, .2, .1, .2...）。

用法：
    python sam_to_fastq.py <input.sam> <output_R1.fq> <output_R2.fq>

示例：
    python sam_to_fastq.py SRR27586274_paired_only.sam fq1.fq fq2.fq
"""

import sys

def sam_to_paired_fastq(sam_path, fq1_path, fq2_path):
    """
    将 SAM 文件转换为两个 FASTQ 文件（R1 和 R2），交替写入。
    注意：此方法假设配对 reads 在文件中按顺序连续出现。
    """
    count = 0        # 总处理行数
    fq1_count = 0    # R1 写入数量
    fq2_count = 0    # R2 写入数量

    try:
        with open(sam_path, 'r') as sam_file, \
             open(fq1_path, 'w') as fq1, \
             open(fq2_path, 'w') as fq2:

            for line in sam_file:
                line = line.strip()
                # 跳过 header 行
                if not line or line.startswith('@'):
                    continue

                parts = line.split('\t', 10)
                if len(parts) < 11:
                    continue  # 字段不足，跳过

                barcode = parts[9][0:16]
                seq = parts[9][20:-1]
                qual = parts[10][20:-1]

                new_name = f"@{parts[0]}:{barcode}"

                # 确保序列和质量值长度一致
                min_len = min(len(seq), len(qual))
                seq = seq[:min_len]
                qual = qual[:min_len]

                # 交替写入 R1 和 R2
                if count % 2 == 0:
                    fq1.write(f"{new_name}\n{seq}\n+\n{qual}\n")
                    fq1_count += 1
                else:
                    fq2.write(f"{new_name}\n{seq}\n+\n{qual}\n")
                    fq2_count += 1

                count += 1

        # 输出统计
        print(f"✅ 处理完成：共 {count} 条 reads")
        print(f"   R1 写入: {fq1_path} （共 {fq1_count} 条）")
        print(f"   R2 写入: {fq2_path} （共 {fq2_count} 条）")

    except FileNotFoundError:
        print(f"❌ 错误：找不到文件 '{sam_path}'。请检查路径是否正确。", file=sys.stderr)
        sys.exit(1)
    except PermissionError:
        print(f"❌ 错误：无法写入输出文件，请检查写入权限。", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"❌ 处理过程中发生未知错误: {e}", file=sys.stderr)
        sys.exit(1)


# ============ 主程序 ============
if __name__ == '__main__':
    # 检查命令行参数数量
    if len(sys.argv) != 4:
        print("❌ 用法错误！", file=sys.stderr)
        print("请按照以下格式运行：", file=sys.stderr)
        print(f"    {sys.argv[0]} <input.sam> <output_R1.fq> <output_R2.fq>", file=sys.stderr)
        print("示例：", file=sys.stderr)
        print(f"    {sys.argv[0]} SRR27586274_paired_only.sam fq1.fq fq2.fq", file=sys.stderr)
        sys.exit(1)

    # 获取命令行参数
    sam_file = sys.argv[1]
    fq1_file = sys.argv[2]
    fq2_file = sys.argv[3]

    # 执行转换
    sam_to_paired_fastq(sam_file, fq1_file, fq2_file)
