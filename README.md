# SV-Insert | 染色体变异种入

A SV Insert project based on Python and High-throughput sequencing. 

一款用Python编写的基于高通测序数据的变异种入程序。

# Introduction | 简介

Write in front: This program is programmed when the first time I learn the Next-Generation-Sequencing and python, and only completed most of the functions. I am a beginner of programming, there is no complete structure design and code separation, and code commenting is not clear. There is nothing call be called algorithm, the most important thing is that the program has no interface, so it is not recommended to use directly after downloading. It is suggested that this procedure should be used when learning about related projects.

写在前面：本程序为初次学习二代测序和Python时编写，只完成了大部分功能，并且本人程序力低下，没有完整的结构设计和代码分离，注释也加的不好，更没有可以被称为算法的东西，最重要的一点是本程序没有接口，所以不建议大家下载之后直接用。建议本程序作为大家对相关项目学习时使用。

The main function of this program is to plant chromosome variations such as tandem repeat variation, inversion variation and translocation variation into the FASTA sequence file. Based on the FASTA file, the fastq file with double terminal read that conforms to the characteristics of the second generation sequencing data is generated according to the setting, and the SAM sequence comparison file is preliminarily screened and analyzed. At the same time, there are some practical small functions, such as BWA command batch generation.

本程序的主要功能为在fasta序列文件中种入串联重复变异、倒位变异、易位变异等染色体变异，以fasta文件为基础按照设定生成生成符合二代测序数据特点的双端read的fastq文件，和初步筛选分析sam序列对比文件。同时附带一些实用的小功能，比如BWA命令批量生成。

Finally, I made up all the English above

最后上面的英文都是我瞎编的
