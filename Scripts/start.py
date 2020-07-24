# coding=utf-8
import os
import my_script
import read_create
import variation_insert
import variant_calling

tandem_repeats_profile_path = r''
inversion_profile_path = r''
translocation_profile_path = r''
default_fasta_path = r''
cmd_file_path = r''
this_path = my_script.create_folder()

while True:
    input_s = ''
    print ('请输入对应数字选择需要进行的程序,输入其他则退出')
    print ('1.变异种入 2.双端read生成 3.变异检测相关程序')
    input_m = raw_input("input:")
    if input_m == '1':
        print ('请输入对应数字选择需要进行的程序,输入其他则返回上层')
        print ('1.串联重复变异 2.倒位变异 3.易位变异')
        input_s = raw_input("input:")
        v_list = []
        while True:
            print ('请输入对应数字选择需要进行的变异编号,输入0则退出输入并开始执行程序')
            input_n = raw_input("input:")
            if input_n == '0':
                break
            else:
                v_list.append(int(input_n))
        for num in v_list:
            if input_s == '1':
                variation_seq_list = variation_insert.tandem_repeats_simulate(tandem_repeats_profile_path, num)
                my_script.fasta_output(variation_seq_list, this_path)
            elif input_s == '2':
                print 2
                variation_seq_list = variation_insert.inversion_simulate(inversion_profile_path, num)
                my_script.fasta_output(variation_seq_list, this_path)
            elif input_s == '3':
                variation_seq_list = variation_insert.translocation_simulate(translocation_profile_path, num)
                my_script.fasta_output(variation_seq_list, this_path)
    elif input_m == '2':
        n_list = []
        while True:
            print ('请输入对应fasta文件的名称,输入0则退出输入并开始执行程序')
            fasta_name = raw_input("input:")
            if fasta_name == '0':
                break
            else:
                n_list.append(fasta_name)
        for fasta_name in n_list:
            if os.path.exists(this_path + '\\' + fasta_name + '.fasta'):
                read_create.reads_base_create(this_path + '\\' + fasta_name + '.fasta')
            elif os.path.exists(default_fasta_path + '\\' + fasta_name + '.fasta'):
                read_create.reads_base_create(default_fasta_path + '\\' + fasta_name + '.fasta')
            else:
                print ('找不到相应文件')
    elif input_m == '3':
        print ('请输入对应数字选择需要进行的程序，输入0则退出输入并开始执行程序')
        print ('1.fastq文件处理命令生成 2.序列对比结果分析')
        input_m = raw_input("input:")
        if input_m == '1':
            read_name = []
            control_name = []
            print ('请输入对应文件夹名称')
            input_path = raw_input("input:")
            print ('请输入对应参考序列名称')
            input_ref = raw_input("input:")
            while True:
                print ('请输入对应读段名称,输入0则退出输入并开始执行程序')
                input_read = raw_input("input:")
                if input_read == '0':
                    break
                else:
                    read_name.append(input_read)
            while True:
                print ('请输入对照组读段名称,输入0则退出输入并开始执行程序')
                input_read = raw_input("input:")
                if input_read == '0':
                    break
                else:
                    control_name.append(input_read)
            cmd_file = open(cmd_file_path + '\\' + input_path + 'cmd_file.txt', 'w+')
            cmd_file.write('cd /media/sf_Py/outdata/' + input_path + '实验集\n')
            cmd_file.write('bwa index -a bwtsw ' + input_path + '.fasta\n')
            variant_calling.cmd_create(input_ref, input_ref, cmd_file, 0)
            for name in read_name:
                variant_calling.cmd_create(name, input_ref, cmd_file, 1)
            for name in control_name:
                variant_calling.cmd_create(name, input_ref, cmd_file, 0)
            cmd_file.close()
        elif input_m == '2':
            pass
    else:
        break
