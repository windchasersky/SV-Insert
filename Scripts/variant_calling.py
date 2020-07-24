# coding=utf-8
import os
import math
insert_k = 2


def cmd_create(seq_name, ref_name, cmd_file, flag):
    cmd_file.write('bwa mem ' + ref_name + '.fasta ' + seq_name + '_read1.fastq ' + seq_name + '_read2.fastq > ' +
                   seq_name + '_pe.sam\n')
    cmd_file.write('samtools view -bS ' + seq_name + '_pe.sam -o ' + seq_name + '_pe.bam\n')
    cmd_file.write('samtools sort ' + seq_name + '_pe.bam  -o ' + seq_name + '_pe.bam.sort\n')
    cmd_file.write('samtools index ' + seq_name + '_pe.bam.sort\n')
    cmd_file.write('samtools mpileup -f ' + ref_name + '.fasta ' + seq_name + '_pe.bam.sort -o ' + seq_name +
                   '_pe.bcf\n')
    cmd_file.write('mkdir ./' + seq_name + '\n')
    if flag == 1:
        cmd_file.write('mv ' + seq_name + '.fasta ./' + seq_name + '\n')
    cmd_file.write('mv ' + seq_name + '_pe.bam ./' + seq_name + '\n')
    cmd_file.write('mv ' + seq_name + '_pe.bam.sort ./' + seq_name + '\n')
    cmd_file.write('mv ' + seq_name + '_pe.bam.sort.bai ./' + seq_name + '\n')
    cmd_file.write('mv ' + seq_name + '_pe.bcf ./' + seq_name + '\n')
    cmd_file.write('mv ' + seq_name + '_pe.sam ./' + seq_name + '\n')
    cmd_file.write('mv ' + seq_name + '_read1.fastq ./' + seq_name + '\n')
    cmd_file.write('mv ' + seq_name + '_read2.fastq ./' + seq_name + '\n')
    cmd_file.write('\n')


# sam遍历筛选
def sam_analysis(sam_path, out_path, length):
    file_obj = open(sam_path)
    file_name = os.path.basename(sam_path).split('.')
    file_name = file_name[0]
    insert_file = open(out_path + '\\' + file_name + '_insert.sam', 'w+')
    flag_file = open(out_path + '\\' + file_name + '_flag.sam', 'w+')
    quality_file = open(out_path + '\\' + file_name + '_quality.sam', 'w+')
    for i in range(2):
        line = file_obj.readline()
        insert_file.write(line)
        flag_file.write(line)
        quality_file.write(line)
    while True:
        line = file_obj.readline()
        if not line:
            break
        line_list = line.split()
        # flag过滤
        flag_list = [0] * 12
        flag = int(line_list[1])
        for i in range(12):
            a = flag % 2
            if a == 1:
                flag_list[i] = 1
            flag = int(flag / 2)
            if flag == 0:
                break
        if flag_list[2] or flag_list[3] or flag_list[8] or flag_list[10] or (flag_list[4] and flag_list[5]):
            flag_file.write(line)
        # 插入长度检测
        if abs(abs(int(line_list[8])) - length) > insert_k * 6 * math.sqrt(length) and flag_list[11] == 0:
            insert_file.write(line)
        quality = int(line_list[4])
        if quality < 30:
            quality_file.write(line)
    file_obj.close()
    insert_file.close()
    flag_file.close()
    quality_file.close()


# 初步分析
def first_analysis(sam_path, insert_file_path, flag_file_path, length, out_path):
    # 初始化
    numbness = 0.01   # 麻木度
    # 文件大小计算
    path_list = [sam_path, insert_file_path, flag_file_path]
    file_size_list = []
    for file_path in path_list:
        count = 1
        file_obj = open(file_path, 'rb')
        while True:
            buffer = file_obj.read(8192 * 1024)
            if not buffer:
                break
            count += buffer.count(b'\n')
        file_size_list.append(int(count)-2)
        file_obj.close()
    sam_file_size = file_size_list[0]
    insert_file_size = file_size_list[1]
    flag_file_size = file_size_list[2]
    # 文件内容提取
    file_obj = open(insert_file_path)
    insert_file_rows = file_obj.readlines()
    del insert_file_rows[0]
    del insert_file_rows[0]
    file_obj.close()

    file_name = os.path.basename(sam_path).split('.')
    first_analysis_file = open(out_path + '\\' + file_name + '_first_analysis.txt', 'w+')
    if insert_file_size > numbness * math.sqrt(sam_file_size) or flag_file_size > numbness * math.sqrt(sam_file_size):
        # 插入长度分析
        insert_flag = 0
        if insert_file_rows:
            insert_length_list = [0, 0]
            for insert_file_row in insert_file_rows:
                length_d = math.fabs(int(insert_file_row.split()[8])) - length
                if length_d > 0:
                    insert_length_list[0] += 1
                else:
                    insert_length_list[1] += 1
            if math.fabs(insert_length_list[0] - insert_length_list[1]) > math.sqrt(insert_file_size):
                if insert_length_list[0] > insert_length_list[1]:
                    insert_flag = 1  # 插入长度过大
                else:
                    insert_flag = 2  # 插入长度过小
            else:
                print '变异过于复杂或插入长度有误'
        # 初步分析
        if math.fabs(insert_file_size - flag_file_size) < numbness * math.sqrt(sam_file_size) / 2:
            first_analysis_file.write('II')
            print '可能包含插入或倒位，进一步分分析请进行局部重对比'
        else:
            if insert_flag == 1:
                first_analysis_file.write('DT')
                print '可能包含缺失或异位变异，进一步分分析请进行局部重对比'
            elif insert_flag == 2:
                first_analysis_file.write('DU')
                print '可能包含重复变异，进一步分分析请进行局部重对比'
            else:
                print '变异过于复杂或插入长度有误'
    else:
        print '未检测出变异'
        return False
    first_analysis_file.close()
