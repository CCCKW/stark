import argparse
import pysam



def reformat_bam(input_sam, output_sam, input_is_bam, output_is_bam, set_bx_flag):

    read_indices = {}
    for align_idx, align in enumerate(input_sam.fetch(until_eof=True)):
        read_id = align.query_name
        read_idx = read_indices.get(read_id, None)
        if read_idx is None:
            read_idx = len(read_indices)
            read_indices[read_id] = read_idx
        if set_bx_flag:
            align.set_tag(tag="BX", value=align.query_name, value_type="Z")
        align.query_name = f"{read_id}:{read_idx}:{align_idx}"
        output_sam.write(align)
    output_sam.close()






def reSAM(input_sam, output_sam, input_is_bam, output_is_bam, set_bx_flag):


    # 打开输入和输出文件
    input_format = "rb" if input_is_bam else "r"
    output_format = "wb" if output_is_bam else "w"

    with pysam.AlignmentFile(input_sam, input_format) as input_sam:
        # 获取输入文件的头部信息，用于初始化输出文件
        header = input_sam.header

        # 打开输出文件，并传递头部信息来初始化它
        with pysam.AlignmentFile(output_sam, output_format, header=header) as output_sam:
            # 调用函数进行重新格式化
            reformat_bam(input_sam, output_sam, input_is_bam, output_is_bam, set_bx_flag)
