# coding=utf-8
# 大量通用代码
import os
import time
import random
import re


# 按照系统时间创建输出文件夹
def create_folder():
    # 引入模块
    localtime = time.strftime('%y%m%d%H%M%S', time.localtime())
    localpath = r'.\output' + localtime
    localpath = localpath.strip()         # 去除首位空格
    localpath = localpath.rstrip("\\")    # 去除尾部 \ 符号

    # 判断路径是否存在
    isExists = os.path.exists(localpath)

    if not isExists:
        # 如果不存在则创建目录,创建目录操作函数
        os.makedirs(localpath)
        print localpath + ' 创建成功'
    else:
        # 如果目录存在则不创建，并提示目录已存在
        print localpath + ' 目录已存在'
    return localpath


# 加载模板链
def load_template_strand(file_path):
    # 初始化
    null_line = 'N' * 80
    seq_list = [0]
    max_len = 1000000

    cut_flag = 1
    loading_flag = 0
    name = ""
    seqlen = 0
    seq = ""
    count = 0
    total_len = 0

    # 获取文件总量
    file_obj = open(file_path, 'rb')
    while True:
        buffer = file_obj.read(8192 * 1024)
        if not buffer:
            break
        count += buffer.count(b'\n')
    total_line = float(count)
    file_obj.close()

    # 加载模板链
    print time.strftime("%m月%d日 %H:%M:%S", time.localtime()) + ' 开始加载模板链'
    template_strand_file = open(file_path, 'r')
    lines = template_strand_file.readlines()
    count = 0
    for line in lines:
        count += 1
        loading_flag = loading_massage(count, total_line, loading_flag)
        line = line.strip('\n')
        if line == "":
            break
        elif line == null_line:
            continue
        elif line[0] == ">":
            name = line.split(' ')[0][1:]
        else:
            seqlen += len(line)
            total_len += len(line)
            seq += line

        if seqlen >= max_len:
            chromosomes = {
                'seqlen': seqlen,
                'seq': seq
            }
            seq_list.append(chromosomes)
            seqlen = 0
            seq = ""
            cut_flag += 1
    chromosomes = {
        'seqlen': seqlen,
        'seq': seq
    }
    seq_list.append(chromosomes)
    seq_list[0] = {
        'name': name,
        'info': cut_flag,
        'seqlen': total_len
    }
    template_strand_file.close()
    print seq_list[0]
    print time.strftime("%m月%d日 %H:%M:%S", time.localtime()) + ' 模板链加载完成'
    return seq_list


# 获取三个数值在中间的
def medium(a, b, c):
    lt = [a, b, c]
    lt.sort()
    return int(lt[1])


# 用标准化信息在指定地点输出fasta文件
def fasta_output(seq_list, path):
    flag = 0
    count = 0
    rest_seq = ''
    isExists = os.path.exists(path)
    if not isExists:
        print time.strftime("%m月%d日 %H:%M:%S", time.localtime()) + ' 未找到相关文件夹'
        return False
    else:
        print time.strftime("%m月%d日 %H:%M:%S", time.localtime()) + ' 开始输出' + seq_list[0]['name'] + '.fasta文件'
        seq_list_len = len(seq_list) -1
        with open(path + '\\' + seq_list[0]['name'] + '.fasta', 'w+') as out_file:
            first_line = '>' + seq_list[0]['name'] + ' ' + seq_list[0]['info'] + \
                         ' seqlen:' + str(seq_list[0]['seqlen']) + '\n'
            out_file.write(first_line)
            del seq_list[0]
            for seq_data in seq_list:
                count += 1
                flag = loading_massage(float(count), seq_list_len, flag)
                seq_data['seq'] = rest_seq + seq_data['seq']
                temp_list = re.findall(r'.{80}', seq_data['seq'])
                rest_seq = seq_data['seq'][len(temp_list) * 80:]
                for temp in temp_list:
                    out_file.write(temp + '\n')
            if rest_seq:
                out_file.write(rest_seq)
        print time.strftime("%m月%d日 %H:%M:%S", time.localtime()) + ' fasta文件输出完成'
        return True


# 用标准信息在指定文件输出fastq信息
def fastq_output(seq_dict, fastq_file):
    fastq_file.write('@' + seq_dict['name'] + '\n')
    fastq_file.write(seq_dict['seq'] + '\n')
    fastq_file.write('+' + '\n')
    fastq_file.write(seq_dict['quality'] + '\n')
    return True


# fasta拆分输出
def fasta_split(fasta_file_path):
    print('读取fasta文件')
    fasta_file = open(fasta_file_path)
    file_usage = 0
    name = 'name'
    seq = []
    for fasta_file_line in fasta_file:
        if fasta_file_line.startswith('>'):
            if not file_usage == 0:
                with open('.\\outdata\\GRCh37\\' + name + '.fasta', 'w+') as out_file:
                    print('正在生成' + str(file_usage))
                    for data in seq:
                        out_file.write(data)
                        seq = []
                    print(str(file_usage) + '生成完成')
            name_line = fasta_file_line.split()
            name = name_line[0].replace('>', '')
            file_usage += 1
        seq.append(fasta_file_line)
    with open('.\\outdata\\GRCh37\\' + name + '.fasta', 'w+') as out_file:
        print(file_usage)
        for data in seq:
            out_file.write(data)
    fasta_file.close()


# 改造fai文件
def mydict_create(faifile_path):
    faifile = open(faifile_path)
    faifile_rows = faifile.readlines()
    faifile.close()
    newfaifile = open(r'.\outdata\My_Homo_sapiens_assembly19.fasta.dict', 'w+')
    contig_top_line = 1
    contig_total_line = 0
    for faifile_row in faifile_rows:
        every_row_cols = faifile_row.split()
        contig_top_line = contig_top_line + contig_total_line
        chrolength = int(every_row_cols[1])
        bp_perline = int(every_row_cols[3])
        contig_total_line = int(chrolength / bp_perline) + 2
        new_row = faifile_row.replace('\n', '\t') + str(contig_total_line) + '\t' + str(contig_top_line) + '\n'
        newfaifile.write(new_row)
    newfaifile.close()


def point_variation(rate, length):
    position = []
    a = 0
    while True:
        a += random.expovariate(rate)
        if a < length and int(round(a)) != 0:
            position.append(int(round(a)))
        else:
            break
    return position


def point_replace(seq, position_list):
    a = ['A', 'T', 'C', 'G']
    for position in position_list:
        r = random.randint(0, 3)
        left = seq[:position - 1]
        right = seq[position:]
        seq = left + a[r] + right
    return seq


def point_insert(seq, position_list):
    a = ['A', 'T', 'C', 'G']
    for position in position_list:
        r = random.randint(0, 3)
        left = seq[:position - 1]
        right = seq[position - 1:]
        seq = left + a[r] + right
    return seq


def point_del(seq, position_list):
    for position in position_list:
        left = seq[:position - 1]
        right = seq[position:]
        seq = left + right
    return seq


def loading_massage(count, total, flag):
    if count/total > 0.25 and flag == 0:
        print time.strftime("%m月%d日 %H:%M:%S", time.localtime()) + ' 完成25%'
        flag += 1
    elif count/total > 0.5 and flag == 1:
        print time.strftime("%m月%d日 %H:%M:%S", time.localtime()) + ' 完成50%'
        flag += 1
    elif count/total > 0.75 and flag == 2:
        print time.strftime("%m月%d日 %H:%M:%S", time.localtime()) + ' 完成75%'
        flag += 1
    return flag

