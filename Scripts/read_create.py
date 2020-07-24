# coding=utf-8
# 读段相关程序
import random
import math
import time
import my_script
from string import maketrans
import os
# 测序参数
insert_length = 300
reads_length = 100
min_reads_length = 50
reads_length_stabilize = 0.8
sequencing_depth = 10
adapter = {
    'name': 'adapter',
    'seq1': 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTTCCGTATCTCGTATGCCGTCTTCTGCTTG',
    'seq2': 'CAAGCAGAAGACGGCATACGAGATACGGAACTGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC',
    'length': 65
}
max_score = 40
pcr_like_switch = False
calling_wrong_switch = True
calling_wrong_rate = 0.03  # 如果是替换就*1.33
read_length_random_switch = True
quality_decrease_switch = True


# 双端read生成
def reads_base_create(seq_path):
    # 初始化读段生成
    pcr = 1
    decrease = 0
    quality_decrease_list = [0] * 20
    flag = 0
    seq_list = my_script.load_template_strand(seq_path)
    fragment_size = sequencing_depth * seq_list[0]['seqlen'] / reads_length / 2

    print time.strftime("%m月%d日 %H:%M:%S", time.localtime()) + ' 为' + seq_list[0]['name'] + '生成' + str(
        fragment_size) + '对read '
    i = 0
    isExists = os.path.exists(seq_path)
    if not isExists:
        print time.strftime("%m月%d日 %H:%M:%S", time.localtime()) + ' 未找到相关文件夹'
        return False
    else:
        read1_file = open(os.path.dirname(seq_path) + '\\' + seq_list[0]['name'] + '_read1.fastq', 'w+')
        read2_file = open(os.path.dirname(seq_path) + '\\' + seq_list[0]['name'] + '_read2.fastq', 'w+')
        # 生成质量递减表
        if quality_decrease_switch:
            quality_decrease_list = quality_decrease(fragment_size, quality_decrease_list)
        # 生成质量列表
        read1_quality_list = quality_score_read1(reads_length, fragment_size)
        read2_quality_list = quality_score_read2(reads_length, fragment_size)
        while True:
            # 生成插入片段
            i += 1
            flag = my_script.loading_massage(i, float(fragment_size), flag)
            cut_dict = insert_seq_create(seq_list, i)
            # 插入片段有N则不记录重新算
            if cut_dict['name'] == 'N':
                i = cut_dict['seqlen']
                continue
            if pcr_like_switch and random.random() < 0.000001:
                pcr = 2
            for p in range(pcr):
                # 生成read
                read1_seq_dict = read1_create(cut_dict)
                read2_seq_dict = read2_create(cut_dict)
                # 生成测序质量
                if quality_decrease_switch:
                    if i > quality_decrease_list[decrease]:
                        decrease += 1
                read1_seq_dict = quality_score_create(read1_seq_dict, decrease, read1_quality_list)
                read2_seq_dict = quality_score_create(read2_seq_dict, decrease, read2_quality_list)
                # 写入fastq文件
                my_script.fastq_output(read1_seq_dict, read1_file)
                my_script.fastq_output(read2_seq_dict, read2_file)
            pcr = 1
            if i >= fragment_size:
                break
        read1_file.close()
        read2_file.close()
    print time.strftime("%m月%d日 %H:%M:%S", time.localtime()) + ' 双端read生成完毕'
    return True


# 插入片段生成
def insert_seq_create(seq_list, i):
    k = 2  # 插入长度波动幅度
    # 插入长度随机与计算
    random_reads_length = int(random.normalvariate(insert_length, k * math.sqrt(insert_length)))
    cut_position1 = random.randint(0, seq_list[0]['seqlen'] - random_reads_length)
    cut_segment = int(cut_position1 / seq_list[1]['seqlen'])
    cut_position1 = cut_position1 - cut_segment * seq_list[1]['seqlen']
    cut_position2 = int(random_reads_length + cut_position1 - seq_list[1]['seqlen'])
    # 生成插入片段
    if cut_position2 > 0:  # 跨越片段时
        cut_seq = seq_list[cut_segment + 1]['seq'][cut_position1:] + \
                  seq_list[cut_segment + 2]['seq'][:cut_position2]
    else:  # 未跨越片段时
        cut_seq = seq_list[cut_segment + 1]['seq'][cut_position1: cut_position1 + random_reads_length]
    # 插入片段有N则重新计算
    if cut_seq.find('N') != -1:
        i -= 1
        cut_dict = {
            'name': 'N',
            'seqlen': i
        }
        return cut_dict
    # 测序错误模拟
    if calling_wrong_switch:
        point = my_script.point_variation(calling_wrong_rate / 3 * 4, random_reads_length)
        if point:
            cut_seq = my_script.point_replace(cut_seq, point)
    # 插入片段输出
    cut_dict = {
        'name': seq_list[0]['name'] + ':' + str(i),
        'seq': cut_seq,
        'seqlen': random_reads_length
    }
    return cut_dict


# read生成
def read1_create(seq_dict):
    this_read_length = reads_length
    # 读长随机
    if read_length_random_switch:
        while True:
            length_ = int(ng_expovariate_random(reads_length, reads_length_stabilize, 0.01))
            if length_ < reads_length - min_reads_length:
                break
        this_read_length = reads_length - length_
    temp_reads_length = min(seq_dict['seqlen'], this_read_length)
    read1 = seq_dict['seq'][:temp_reads_length]
    # 插入片段过小时加入接头序列补齐
    if seq_dict['seqlen'] < this_read_length:
        compensate_length = this_read_length - seq_dict['seqlen']
        read1 = read1 + adapter['seq1'][: compensate_length]
    # 信息记录
    read1_dict = {
        'name': seq_dict['name'] + '/1',
        'seq': read1,
        'length': this_read_length
    }
    return read1_dict


def read2_create(seq_dict):
    this_read_length = reads_length
    in_tab = "ATCG"
    out_tab = "TAGC"
    tran_tab = maketrans(in_tab, out_tab)
    # 读长随机
    if read_length_random_switch:
        while True:
            length_ = int(ng_expovariate_random(reads_length, reads_length_stabilize, 0.01))
            if length_ < reads_length - min_reads_length:
                break
        this_read_length = reads_length - length_
    temp_reads_length = min(seq_dict['seqlen'], this_read_length)
    read2 = seq_dict['seq'][:-temp_reads_length - 1:-1]
    read2 = read2.translate(tran_tab)
    if seq_dict['seqlen'] < this_read_length:
        compensate_length = this_read_length - seq_dict['seqlen']
        read2 = read2 + adapter['seq2'][: compensate_length]
    # 信息记录
    read2_dict = {
        'name': seq_dict['name'] + '/2',
        'seq': read2,
        'length': this_read_length
    }
    return read2_dict


# 测序质量递减表生成
def quality_decrease(size, decrease_list):
    a = 0
    for i in range(10000):
        while True:
            decrease = int(ng_expovariate_random(20, 0.9, 0.1))
            if decrease < 20:
                break
        decrease_list[decrease] += 1
    for i in range(20):
        a += decrease_list[i] / 10000.0 * size
        decrease_list[i] = int(a)
    return decrease_list


# 质量分数生成
def quality_score_create(read_seq_dict, decrease, quality_list):
    # 生成质量
    score_string = ''
    random_list = random.randint(0, len(quality_list) - 1)
    for i in range(read_seq_dict['length'] - 1):
        score_string += chr(quality_list[random_list][i] - decrease)
    score_string += chr(quality_list[random_list][read_seq_dict['length'] - 1] - 3)
    read_dict = {
        'name': read_seq_dict['name'],
        'seq': read_seq_dict['seq'],
        'quality': score_string
    }
    return read_dict


# 质量列表生成
def quality_score_read1(length, size):
    c = -0.084117647
    k = 0.15
    b = 34 - c * 99
    quality_list = []
    for i in range(int(math.sqrt(size))):
        quality = []
        for position in range(length):
            if position < 3:
                score = 33
            elif position < 14:
                score = random.normalvariate(36.6 - 0.01 * position, k * math.sqrt(position))
            else:
                score = random.normalvariate(c * position + b, k * math.sqrt(position))
            score = min(max_score, score)
            quality_score = int(score) + 33
            quality.append(quality_score)
        quality_list.append(quality)
    return quality_list


def quality_score_read2(length, size):
    c = -0.0705882352941176
    k = 0.3
    b = 33.6 - c * 99
    quality_list = []
    for i in range(int(math.sqrt(size))):
        quality = []
        for position in range(length):
            if position < 3:
                score = 32
            elif position < 14:
                score = random.normalvariate(35.62 - 0.01 * position, k * math.sqrt(position))
            else:
                score = random.normalvariate(c * position + b, k * math.sqrt(position))
            score = min(max_score, score)
            quality_score = int(score) + 33
            quality.append(quality_score)
        quality_list.append(quality)
    return quality_list


# 叠加指数分布
def ng_expovariate_random(lam, stabilize, k):
    a = random.expovariate(1 / float(lam))
    b = random.random()
    if b < stabilize:
        a = random.expovariate(k * lam)
    return a
