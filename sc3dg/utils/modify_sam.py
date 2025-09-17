#!/usr/bin/env python3

def modify_sam_file(input_sam, output_sam):
    """
    修改 SAM 文件：
    - 将 QNAME 中的 'barcode:CB:UMI' 拆分
    - 保留 barcode 在 QNAME
    - 将 CB 和 UMI 附加到 SEQ 和 QUAL 字段，用 '@**@' 分隔

    示例：
        原始 QNAME: ABC:CTGAACGGACTC:ATGCCTAA
        修改后:
            QNAME: ABC
            SEQ:   original_seq@**@CTGAACGGACTC:ATGCCTAA
            QUAL:  original_qual@**@CTGAACGGACTC:ATGCCTAA
    """
    with open(input_sam, 'r') as fin, open(output_sam, 'w') as fout:
        for line_num, line in enumerate(fin, 1):
            line = line.rstrip('\n')  # 保留原始换行处理

            # 如果是 header 行（以 @ 开头），直接输出
            if line.startswith('@'):
                fout.write(line + '\n')
                continue

            # 跳过空行
            if not line.strip():
                continue

            components = line.split('\t')

            # 确保至少有 11 列（SAM 要求至少 11 列）
            if len(components) < 11:
                print(f"⚠️ 第 {line_num} 行字段不足，跳过: {line}", file=sys.stderr)
                continue

            qname = components[0]
            seq = components[9]
            qual = components[10]

            # 拆分 QNAME: 格式应为 barcode:CB:UMI
            subparts = qname.split(':', 2)  # 最多拆成3部分
            if len(subparts) != 3:
                print(f"⚠️ 第 {line_num} 行 QNAME 格式错误（期望 'bc:CB:UMI'）: {qname}", file=sys.stderr)
                fout.write(line + '\n')  # 保留原行
                continue

            # 提取各部分
            barcode = subparts[0]  # 如 SRR27586274
            readseq = subparts[1]  # 如 CTGAACGGACTC
            readqua = subparts[2]  # 如 ATGCCTAA

            # 修改 components
            components[0] = barcode
            components[9] = seq + "@**@" + readseq
            components[10] = qual + "@**@" + readqua
            components[5] = str(len(components[9]))+'M'
            # 写入输出文件（带换行）
            fout.write('\t'.join(components) + '\n')

    print(f"✅ 处理完成！输出文件: {output_sam}")


# ============ 主程序 ============
if __name__ == '__main__':
    import sys

    if len(sys.argv) != 3:
        print("📌 用法: python modify_sam.py <input.sam> <output.sam>", file=sys.stderr)
        print("示例: python modify_sam.py SRR27586274_R1_BC.sam SRR27586274_R1_BC_modified.sam", file=sys.stderr)
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    try:
        modify_sam_file(input_file, output_file)
    except FileNotFoundError:
        print(f"❌ 错误：找不到输入文件 '{input_file}'", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"❌ 处理过程中发生错误: {e}", file=sys.stderr)
        sys.exit(1)