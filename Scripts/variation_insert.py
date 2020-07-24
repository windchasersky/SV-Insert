# coding=utf-8
# 变异种入相关代码
# 串联变异配置文件 1：变异代号 2：变异染色体名称 3：变异位置 4：重复段长度
# 倒位变异配置文件 1：变异代号 2：变异染色体名称 3：变异位置 4：倒位段长度
# 易位变异配置文件 1：变异代号 2：变异染色体名称 3：变异位置1 4：变异位置2 5：易位段长度 6：易位方向（0前面插入后面，1后面插入前面）
import math
import os
import random
import time
import my_script


# 通过变异配置文件为模板链种入串联重复变异
def tandem_repeats_simulate(tandem_repeats_profile_path, num):
    # 初始化
    repeat_time_min = 3
    repeat_time_max = 10
    random_switch = True
    # 载入配置信息
    print time.strftime("%m月%d日 %H:%M:%S", time.localtime()) + ' 加载变异配置信息'
    variation_profile = open(tandem_repeats_profile_path)
    variation_profile_rows = variation_profile.readlines()
    variation_profile_row = variation_profile_rows[num]
    variation_profile.close()
    variation_profile_data = variation_profile_row.split()
    tandem_repeats_list = {
        'vid': variation_profile_data[0],
        'name': variation_profile_data[1],
        'position': int(variation_profile_data[2]),
        'length': int(variation_profile_data[3]),
        'times': random.randint(repeat_time_min, repeat_time_max)
    }
    # 种入变异
    print time.strftime("%m月%d日 %H:%M:%S", time.localtime()) + ' 开始' + tandem_repeats_list['vid'] + '变异'
    template_strand_file_path = '.\\outdata\\GRCh37\\' + tandem_repeats_list['name'] + '.fasta'
    is_exists = os.path.exists(template_strand_file_path)
    template_strand = my_script.load_template_strand(template_strand_file_path)  # 载入模板链
    if tandem_repeats_list['position'] >= template_strand[0]['seqlen']:
        print time.strftime("%m月%d日 %H:%M:%S", time.localtime()) + ' ' + tandem_repeats_list['vid'] + '变异位点信息有误'
        return False
    if is_exists:
        # 随机化变异参数
        if random_switch:
            random_repeat_length = int(
                random.normalvariate(tandem_repeats_list['length'], math.sqrt(tandem_repeats_list['length'])))
            tandem_repeats_list['length'] = my_script.medium(50, random_repeat_length, template_strand[0]['seqlen'])
            random_position = int(
                random.normalvariate(tandem_repeats_list['position'], math.sqrt(template_strand[0]['seqlen'])))
            tandem_repeats_list['position'] = my_script.medium(tandem_repeats_list['length'], random_position,
                                                               template_strand[0]['seqlen'])
        # 计算变异参数
        v_segment1 = int(tandem_repeats_list['position'] / template_strand[1]['seqlen'])
        tandem_repeats_list['position'] = tandem_repeats_list['position'] - v_segment1 * template_strand[1]['seqlen']
        v_position2 = tandem_repeats_list['length'] - tandem_repeats_list['position']
        # 生成变异链
        print time.strftime("%m月%d日 %H:%M:%S", time.localtime()) + ' 开始' + tandem_repeats_list['vid'] + '变异种入'
        # 生成重复片段
        if v_position2 > 0:  # 跨越片段时
            v_segment2 = int(v_position2 / template_strand[1]['seqlen'])
            v_position2 = v_segment2 * template_strand[1]['seqlen'] - v_position2
            repeat_seq = template_strand[v_segment1 - v_segment2]['seq'][v_position2:]
            if v_segment2 >= 1:
                for i in range(v_segment2 - 1):
                    repeat_seq += template_strand[v_segment1 - v_segment2 + 1 + i]['seq']
            repeat_seq += template_strand[v_segment1 + 1]['seq'][:tandem_repeats_list['position']]
        else:  # 未跨越片段时
            repeat_seq = template_strand[v_segment1 + 1]['seq'][
                         tandem_repeats_list['position'] - tandem_repeats_list['length']:
                         tandem_repeats_list['position']]
        # 生成变异片段
        seq_medium = ''
        for i in range(tandem_repeats_list['times']):
            seq_medium += repeat_seq
        seq_left = template_strand[v_segment1 + 1]['seq'][:tandem_repeats_list['position']]
        seq_right = template_strand[v_segment1 + 1]['seq'][tandem_repeats_list['position']:]
        variation_seq = seq_left + seq_medium + seq_right
        # 更新列表信息
        variation_seq_dict = {
            'id': v_segment1 + 1,
            'seqlen': int(template_strand[v_segment1 + 1]['seqlen']) + (tandem_repeats_list['length'] *
                                                                        tandem_repeats_list['times']),
            'seq': variation_seq
        }
        template_strand[v_segment1 + 1] = variation_seq_dict
        variation_message_dict = {
            'name': tandem_repeats_list['vid'],
            'info': 'chroname:' + tandem_repeats_list['name'] + ' ' +
                    'vposition:' + str(tandem_repeats_list['position']) + ' ' +
                    'vlength:' + str(tandem_repeats_list['length']) + ' ' +
                    'vtimes:' + str(tandem_repeats_list['times']),
            'seqlen': int(template_strand[0]['seqlen']) + tandem_repeats_list['length'] * tandem_repeats_list['times'],
        }
        template_strand[0] = variation_message_dict
        print time.strftime("%m月%d日 %H:%M:%S", time.localtime()) + ' ' + tandem_repeats_list['vid'] + '变异完成'
        return template_strand
    else:
        print time.strftime("%m月%d日 %H:%M:%S", time.localtime()) + ' ' + tandem_repeats_list['vid'] + '变异未匹配到相关文件'
        return False


# 通过变异配置文件为模板链种入倒位变异
def inversion_simulate(inversion_profile_path, num):
    # 初始化
    random_switch = True
    # 载入配置信息
    print time.strftime("%m月%d日 %H:%M:%S", time.localtime()) + ' 加载变异配置信息'
    variation_profile = open(inversion_profile_path)
    variation_profile_rows = variation_profile.readlines()
    variation_profile_row = variation_profile_rows[num]
    variation_profile.close()
    variation_profile_data = variation_profile_row.split()
    inversion_list = {
        'vid': variation_profile_data[0],
        'name': variation_profile_data[1],
        'position': int(variation_profile_data[2]),
        'length': int(variation_profile_data[3]),
    }
    # 种入变异
    print time.strftime("%m月%d日 %H:%M:%S", time.localtime()) + ' 开始' + inversion_list['vid'] + '变异'
    template_strand_file_path = '.\\outdata\\GRCh37\\' + inversion_list['name'] + '.fasta'
    is_exists = os.path.exists(template_strand_file_path)
    template_strand = my_script.load_template_strand(template_strand_file_path)  # 载入模板链
    if inversion_list['position'] >= template_strand[0]['seqlen']:
        print time.strftime("%m月%d日 %H:%M:%S", time.localtime()) + ' ' + inversion_list['vid'] + '变异位点信息有误'
        return False
    if is_exists:
        # 随机化变异参数
        if random_switch:
            random_repeat_length = int(
                random.normalvariate(inversion_list['length'], math.sqrt(inversion_list['length'])))
            inversion_list['length'] = my_script.medium(50, random_repeat_length, template_strand[0]['seqlen'])
            random_position = int(
                random.normalvariate(inversion_list['position'], math.sqrt(template_strand[0]['seqlen'])))
            inversion_list['position'] = my_script.medium(inversion_list['length'], random_position,
                                                          template_strand[0]['seqlen'])
        # 计算变异参数
        v_segment1 = int(inversion_list['position'] / template_strand[1]['seqlen'])
        inversion_list['position'] = inversion_list['position'] - v_segment1 * template_strand[1]['seqlen']
        v_position2 = inversion_list['length'] - inversion_list['position']
        # 生成变异链
        print time.strftime("%m月%d日 %H:%M:%S", time.localtime()) + ' 开始' + inversion_list['vid'] + '变异种入'
        if v_position2 > 0:  # 跨越片段时
            v_segment2 = int(v_position2 / template_strand[1]['seqlen'])
            v_position2 = v_segment2 * template_strand[1]['seqlen'] - v_position2
            # 前
            inversion_seq = template_strand[v_segment1 - v_segment2]['seq'][v_position2:]
            inversion_seq = inversion_seq[::-1]
            seq_left = template_strand[v_segment1 - v_segment2]['seq'][:v_position2]
            variation_seq = seq_left + inversion_seq
            template_strand[v_segment1 - v_segment2]['seq'] = variation_seq
            # 中
            if v_segment2 >= 1:
                for i in range(v_segment2 - 1):
                    inversion_seq = template_strand[v_segment1 - v_segment2 + 1 + i]['seq'][::-1]
                    template_strand[v_segment1 - v_segment2 + 1 + i]['seq'] = inversion_seq
            # 后
            inversion_seq = template_strand[v_segment1 + 1]['seq'][:inversion_list['position']]
            inversion_seq = inversion_seq[::-1]
            seq_right = template_strand[v_segment1 + 1]['seq'][inversion_list['position']:]
            variation_seq = inversion_seq + seq_right
            template_strand[v_segment1 + 1]['seq'] = variation_seq
        else:  # 未跨越片段时
            inversion_seq = template_strand[v_segment1 + 1]['seq'][
                            inversion_list['position'] - inversion_list['length']:inversion_list['position']]
            inversion_seq = inversion_seq[::-1]
            seq_left = template_strand[v_segment1 + 1]['seq'][:inversion_list['position'] - inversion_list['length']]
            seq_right = template_strand[v_segment1 + 1]['seq'][inversion_list['position']:]
            variation_seq = seq_left + inversion_seq + seq_right
            template_strand[v_segment1 + 1]['seq'] = variation_seq
        # 更新列表信息
        v_dict = {
            'name': inversion_list['vid'],
            'info': 'chroname:' + inversion_list['name'] + ' ' +
                    'vposition:' + str(inversion_list['position'] + v_segment1 * template_strand[1]['seqlen']) + ' ' +
                    'vlength:' + str(inversion_list['length']),
            'seqlen': int(template_strand[0]['seqlen'])
        }
        template_strand[0] = v_dict
        print time.strftime("%m月%d日 %H:%M:%S", time.localtime()) + ' ' + inversion_list['vid'] + '变异完成'
        return template_strand
    else:
        print time.strftime("%m月%d日 %H:%M:%S", time.localtime()) + ' ' + inversion_list['vid'] + '变异未匹配到相关文件'
        return False


# 通过变异配置文件为模板链种入易位变异
def translocation_simulate(translocation_profile_path, num):
    # 初始化
    min_distance = 50
    random_switch = True
    # 载入配置信息
    print time.strftime("%m月%d日 %H:%M:%S", time.localtime()) + ' 加载变异配置信息'
    variation_profile = open(translocation_profile_path)
    variation_profile_rows = variation_profile.readlines()
    variation_profile_row = variation_profile_rows[num]
    variation_profile.close()
    variation_profile_data = variation_profile_row.split()
    translocation_list = {
        'vid': variation_profile_data[0],
        'name': variation_profile_data[1],
        'position1': int(variation_profile_data[2]),
        'position2': int(variation_profile_data[3]),
        'length': int(variation_profile_data[4]),
        'direction': int(variation_profile_data[5])
    }

    # 种入变异
    print time.strftime("%m月%d日 %H:%M:%S", time.localtime()) + ' 开始' + translocation_list['vid'] + '变异'
    template_strand_file_path = '.\\outdata\\GRCh37\\' + translocation_list['name'] + '.fasta'
    is_exists = os.path.exists(template_strand_file_path)
    template_strand = my_script.load_template_strand(template_strand_file_path)  # 载入模板链
    if translocation_list['position1'] >= template_strand[0]['seqlen'] or \
            translocation_list['position2'] >= template_strand[0]['seqlen']:
        print time.strftime("%m月%d日 %H:%M:%S", time.localtime()) + ' ' + translocation_list['vid'] + '变异位点信息有误'
        return False
    if is_exists:
        position1 = min(translocation_list['position1'], translocation_list['position2'])
        position2 = max(translocation_list['position1'], translocation_list['position2'])
        # 随机化变异参数
        if random_switch:
            while True:
                random_position1 = int(
                    random.normalvariate(translocation_list['position1'], math.sqrt(template_strand[0]['seqlen'])))
                random_position1 = my_script.medium(0, random_position1, template_strand[0]['seqlen'])
                random_position2 = int(
                    random.normalvariate(translocation_list['position2'], math.sqrt(template_strand[0]['seqlen'])))
                random_position2 = my_script.medium(0, random_position2, template_strand[0]['seqlen'])
                position1 = min(random_position1, random_position2)
                position2 = max(random_position1, random_position2)
                if position2 - position1 >= min_distance:
                    break
        # 计算变异参数
        translocation_list['position1'] = position1
        translocation_list['position2'] = position2
        distance = translocation_list['position2'] - translocation_list['position1']
        random_repeat_length = int(
            random.normalvariate(translocation_list['length'], math.sqrt(translocation_list['length'])))
        translocation_list['length'] = my_script.medium(min_distance, random_repeat_length, distance)
        v_segment1 = int(translocation_list['position1'] / template_strand[1]['seqlen'])
        v_segment2 = int(translocation_list['position2'] / template_strand[1]['seqlen'])
        translocation_list['position1'] = translocation_list['position1'] - v_segment1 * template_strand[1]['seqlen']
        translocation_list['position2'] = translocation_list['position2'] - v_segment2 * template_strand[1]['seqlen']
        if translocation_list['direction'] == 0:
            v_position2 = translocation_list['length'] + translocation_list['position1'] - template_strand[1]['seqlen']
            v_segment_2 = int(v_position2 / template_strand[1]['seqlen'])
        elif translocation_list['direction'] == 1:
            v_position2 = translocation_list['length'] - translocation_list['position2']
            v_segment_2 = int(v_position2 / template_strand[1]['seqlen'])
        else:
            print time.strftime("%m月%d日 %H:%M:%S", time.localtime()) + ' ' + translocation_list['vid'] + '易位方向错误'
            return False
        # 生成变异链并更新列表信息
        print time.strftime("%m月%d日 %H:%M:%S", time.localtime()) + ' 开始' + translocation_list['vid'] + '变异种入'
        if translocation_list['direction'] == 0:   # 前面的片段易位到后面
            if v_position2 > 0:  # 跨越片段时
                v_position2 = v_position2 - v_segment_2 * template_strand[1]['seqlen']
                # 前
                translocation_seq = template_strand[v_segment1 + 1]['seq'][translocation_list['position1']:]
                template_strand[v_segment1 + 1]['seq'] = template_strand[v_segment1 + 1]['seq'][
                                                         :translocation_list['position1']]
                template_strand[v_segment1 + 1]['seqlen'] = len(template_strand[v_segment1 + 1]['seq'])
                # 中
                if v_segment_2 >= 1:
                    for i in range(v_segment_2 - 1):
                        translocation_seq += template_strand[v_segment1 + i + 2]['seq']
                        del template_strand[v_segment1 + i + 2]
                # 后
                translocation_seq += template_strand[v_segment1 + 2]['seq'][:v_position2]
                # 目标位点插入易位片段
                left_seq = template_strand[v_segment2 + 1 - v_segment_2]['seq'][:translocation_list['position2']]
                right_seq = template_strand[v_segment2 + 1 - v_segment_2]['seq'][translocation_list['position2']:]
                template_strand[v_segment2 + 1 - v_segment_2]['seq'] = left_seq + translocation_seq + right_seq
                template_strand[v_segment2 + 1 - v_segment_2]['seqlen'] = len(template_strand[v_segment2 + 1 -
                                                                                              v_segment_2]['seq'])
                template_strand[v_segment1 + 2]['seq'] = template_strand[v_segment1 + 2]['seq'][v_position2:]
                template_strand[v_segment1 + 2]['seqlen'] = len(template_strand[v_segment1 + 2]['seq'])
            else:  # 未跨越片段时
                translocation_seq = template_strand[v_segment1 + 1]['seq'][translocation_list['position1']:
                                                                           translocation_list['position1'] +
                                                                           translocation_list['length']]
                left_seq = template_strand[v_segment2 + 1]['seq'][:translocation_list['position2']]
                right_seq = template_strand[v_segment2 + 1]['seq'][translocation_list['position2']:]
                template_strand[v_segment2 + 1]['seq'] = left_seq + translocation_seq + right_seq
                template_strand[v_segment2 + 1]['seqlen'] = len(template_strand[v_segment2 + 1]['seq'])
                left_seq = template_strand[v_segment1 + 1]['seq'][:translocation_list['position1']]
                right_seq = template_strand[v_segment1 + 1]['seq'][translocation_list['position1'] +
                                                                   translocation_list['length']:]
                template_strand[v_segment1 + 1]['seq'] = left_seq + right_seq
                template_strand[v_segment1 + 1]['seqlen'] = len(template_strand[v_segment1 + 1]['seq'])
        elif translocation_list['direction'] == 1:   # 后面的片段易位到前面
            if v_position2 > 0:  # 跨越片段时
                v_position2 = v_position2 - v_segment_2 * template_strand[1]['seqlen']
                # 前
                translocation_seq = template_strand[v_segment2 - v_segment_2]['seq'][v_position2:]
                template_strand[v_segment2 - v_segment_2]['seq'] = template_strand[v_segment2 + 1 - v_segment_2][
                                                                           'seq'][:v_position2]
                template_strand[v_segment2 - v_segment_2]['seqlen'] = len(
                    template_strand[v_segment2 + 1 - v_segment_2]['seq'])
                # 中
                if v_segment_2 >= 1:
                    for i in range(v_segment_2 - 1):
                        translocation_seq += template_strand[v_segment2 - v_segment_2 + 1]['seq']
                        del template_strand[v_segment2 - v_segment_2 + 1]
                # 后
                translocation_seq += template_strand[v_segment2 + 1 - v_segment_2]['seq'][
                                     :translocation_list['position2']]
                template_strand[v_segment2 + 1 - v_segment_2]['seq'] = template_strand[v_segment1 + 1 - v_segment_2][
                                                                           'seq'][translocation_list['position2']:]
                template_strand[v_segment2 + 1 - v_segment_2]['seqlen'] = len(template_strand[v_segment1 + 1 -
                                                                                              v_segment_2]['seq'])
                # 目标位点插入易位片段
                left_seq = template_strand[v_segment1 + 1]['seq'][:translocation_list['position1']]
                right_seq = template_strand[v_segment1 + 1]['seq'][translocation_list['position1']:]
                template_strand[v_segment1 + 1]['seq'] = left_seq + translocation_seq + right_seq
                template_strand[v_segment1 + 1]['seqlen'] = len(template_strand[v_segment2 + 1 - v_segment_2]['seq'])
            else:  # 未跨越片段时
                translocation_seq = template_strand[v_segment2 + 1]['seq'][translocation_list['position2'] -
                                                                           translocation_list['length']:
                                                                           translocation_list['position2']]
                left_seq = template_strand[v_segment2 + 1]['seq'][
                           :translocation_list['position2'] - translocation_list['length']]
                right_seq = template_strand[v_segment2 + 1]['seq'][translocation_list['position2']:]
                template_strand[v_segment2 + 1]['seq'] = left_seq + right_seq
                template_strand[v_segment2 + 1]['seqlen'] = len(template_strand[v_segment1 + 1]['seq'])
                left_seq = template_strand[v_segment1 + 1]['seq'][:translocation_list['position1']]
                right_seq = template_strand[v_segment1 + 1]['seq'][translocation_list['position1']:]
                template_strand[v_segment1 + 1]['seq'] = left_seq + translocation_seq + right_seq
                template_strand[v_segment1 + 1]['seqlen'] = len(template_strand[v_segment2 + 1]['seq'])

        v_dict = {
            'name': translocation_list['vid'],
            'info': 'chroname:' + translocation_list['name'] + ' ' +
                    'vposition1:' + str(translocation_list['position1']) + ' ' +
                    'vposition2:' + str(translocation_list['position2']) + ' ' +
                    'vlength:' + str(translocation_list['length']) + ' ' +
                    'vdirection:' + str(translocation_list['direction']),
            'seqlen': int(template_strand[0]['seqlen'])
        }
        template_strand[0] = v_dict
        print time.strftime("%m月%d日 %H:%M:%S", time.localtime()) + ' ' + translocation_list['vid'] + '变异完成'
        return template_strand
    else:
        print time.strftime("%m月%d日 %H:%M:%S", time.localtime()) + ' ' + translocation_list['vid'] + '变异未匹配到相关文件'
        return False
