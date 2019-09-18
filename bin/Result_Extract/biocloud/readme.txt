+------------------------------------------------------------------------------+
|                    转录组项目信息分析结果说明文档                            |
+------------------------------------------------------------------------------+

目录结构及文件说明：
********************************************************************************
BMKXXXXXX-XXX_Transcriptome_final/
|-- cleandata        #测序数据评估目录
|   |-- AllSample_GC_Q.stat  #测序数据评估统计表
|   `-- PNG                  #以Cycle为单位对测序数据进行评估作图的目录
|       |-- Sample1.acgtn.png          #样品1测序数据碱基含量分布图
|       |-- Sample1.cycleQ.png         #样品1测序数据平均质量值分布图
|       |-- Sample1.quality.png        #样品1测序数据碱基测序错误分布图
|       |-- Sample1.rawDataStat.png    #样品1测序数据数分布图
|       `-- ... ...
|-- Assembly       #转录组组装结果及统计目录
|   |-- Assemby1             #第1套转录组组装结果目录
|   |   |-- contigs                    #Contig组装结果目录
|   |   |   |-- Assemby1.contigs.distribution.png          #Contig序列长度分布图（png格式）
|   |   |   |-- Assemby1.contigs.distribution.svg          #Contig序列长度分布图（svg格式）
|   |   |   |-- Assemby1.contigs.fa                        #Contig序列文件（FASTA格式）
|   |   |   |-- Assemby1.contigs.n50                       #Contig组装长度统计文件，如N50等
|   |   |   |-- Assemby1.contigs.stat.info.xls             #Contig序列长度分布统计表
|   |   |   `-- Assemby1.contigs.stat.xls                  #Contig组装结果统计表
|   |   |-- Transcripts                #转录本组装结果目录
|   |   |   |-- Assemby1.Transcripts.distribution.png      #转录本序列长度分布图（png格式）
|   |   |   |-- Assemby1.Transcripts.distribution.svg      #转录本序列长度分布图（svg格式）
|   |   |   |-- Assemby1.Transcripts.fa                    #转录本序列文件（FASTA格式）
|   |   |   |-- Assemby1.Transcripts.n50                   #转录本组装长度统计文件，如N50等
|   |   |   |-- Assemby1.Transcripts.stat.info.xls         #转录本序列长度分布统计表
|   |   |   `-- Assemby1.Transcripts.stat.xls              #转录本组装结果统计表
|   |   `-- Unigenes                   #Unigene序列提取结果目录
|   |       |-- Assemby1.Unigene.Cluster.stat.xls          #基因可变剪接统计表
|   |       |-- Assemby1.Unigene_Trans_ID.list             #Unigene和转录本编号对应表
|   |       |-- Assemby1.Unigenes.Cluster.xls              #Unigene及其对应的转录本集文件
|   |       |-- Assemby1.Unigenes.distribution.png         #Unigene序列长度分布图（png格式）
|   |       |-- Assemby1.Unigenes.distribution.svg         #Unigene序列长度分布图（svg格式）
|   |       |-- Assemby1.Unigenes.fa                       #Unigene序列文件（FASTA格式）
|   |       |-- Assemby1.Unigenes.n50                      #Unigene长度统计文件，如N50等
|   |       |-- Assemby1.Unigenes.stat.info.xls            #Unigene长度分布统计表
|   |       `-- Assemby1.Unigenes.stat.xls                 #Unigene提取结果统计表
|   `-- ... ...              #（非合并组装的项目会有1套以上的转录组组装结果）
|-- Unigene        #构建的Unigene库及其功能注释、SSR、CDS预测分析结果目录
|   |-- SpeciesName.Unigene.distribution.png     #目的物种的Unigene库中序列长度分布图（png格式）
|   |-- SpeciesName.Unigene.fa                   #目的物种的Unigene库文件（FASTA格式）
|   |-- SpeciesName.Unigene.stat.info.xls        #目的物种的Unigene库中序列长度分布统计表
|   |-- SpeciesName.Unigene.stat.xls             #目的物种的Unigene库中序列统计表
|   |-- Unigene_Anno         #Unigene功能注释分析结果目录
|   |   |-- All_Database_annotation.xls                    #功能注释、富集和统计整合文件（excel格式）
|   |   |-- Function_Annotation.stat.xls                   #注释数目统计表
|   |   |-- Integrated_Function.annotation.xls             #功能注释整合文件
|   |   |-- Kegg_map                                       #注释到的KEGG通路图目录
|   |   |   |-- ko00010.png
|   |   |   |-- ... ...
|   |   |   `-- ko04626.png
|   |   |-- SpeciesName.Unigene.fa.Cog.cluster.png         #COG注释分类统计图
|   |   |-- SpeciesName.Unigene.fa.Cog.cluster.stat        #COG注释分类统计表
|   |   |-- SpeciesName.Unigene.fa.Cog_class.txt           #COG数据库注释结果文件
|   |   |-- SpeciesName.Unigene.fa.GO.anno.txt             #GO数据库注释结果文件
|   |   |-- SpeciesName.Unigene.fa.GO.list.txt             #Unigene编号与注释到的GO节点编号对应表
|   |   |-- SpeciesName.Unigene.fa.GO.png                  #GO二级节点注释统计图
|   |   |-- SpeciesName.Unigene.fa.GO_enrichment.stat.xls  #GO二级节点注释统计表
|   |   |-- SpeciesName.Unigene.fa.GO_tree.stat.xls        #GO节点与注释到该节点的Unigene编号对应表
|   |   |-- SpeciesName.Unigene.fa.Kegg.globe              #KEGG全局通路与注释到该全局通路的Unigene编号对应表
|   |   |-- SpeciesName.Unigene.fa.Kegg.ko                 #KEGG数据库注释结果文件
|   |   |-- SpeciesName.Unigene.fa.Kegg.pathway            #KEGG通路与注释到该通路的Unigene编号对应表
|   |   |-- SpeciesName.Unigene.fa.Kog.cluster.png         #KOG注释分类统计图
|   |   |-- SpeciesName.Unigene.fa.Kog.cluster.stat        #KOG注释分类统计表
|   |   |-- SpeciesName.Unigene.fa.Kog_class.txt           #KOG数据库注释结果文件
|   |   |-- SpeciesName.Unigene.fa.Pfam.anno.txt           #Pfam数据库注释结果文件
|   |   |-- SpeciesName.Unigene.fa.Swissprot.anno.txt      #Swiss-Prot数据库注释结果文件
|   |   |-- SpeciesName.Unigene.fa.nr.anno.txt             #nr数据库注释结果文件
|   |   |-- SpeciesName.Unigene.fa.nr.lib.png              #nr数据库中注释到的物种分布统计饼图
|   |   `-- SpeciesName.Unigene.fa.nr.lib.stat             #nr数据库中注释到的物种分布统计表
|   |-- Unigene_CDS          #Unigene的CDS预测分析结果目录
|   |   |-- SpeciesName.Unigene.cds.distribution.png       #编码区序列长度分布图（png格式）
|   |   |-- SpeciesName.Unigene.cds.distribution.svg       #编码区序列长度分布图（svg格式）
|   |   |-- SpeciesName.Unigene.cds.fa                     #编码区序列文件（FASTA格式）
|   |   |-- SpeciesName.Unigene.cds.stat.xls               #编码区序列长度分布统计表
|   |   `-- SpeciesName.Unigene.pep.fa                     #氨基酸序列文件（FASTA格式）
|   `-- Unigene_SSR          #1Kb以上Unigene的SSR预测及其引物设计分析结果目录
|       |-- SpeciesName.Unigene.1000.fa                    #1Kb以上的Unigene序列文件
|       |-- SpeciesName.Unigene.1000.fa.ssr.density.png    #SSR密度分布图
|       |-- SpeciesName.Unigene.1000.fa.ssr.density.xls    #SSR密度分布统计表
|       |-- SpeciesName.Unigene.1000.fa.ssr.p3out          #SSR引物设计软件Primer3输出文件
|       |-- SpeciesName.Unigene.1000.fa.stat.xls           #各类型SSR数目统计表
|       |-- SpeciesName.Unigene.1000.fa.SSR.result.xls     #SSR结果和引物设计结果
|       `-- SpeciesName.Unigene.1000.fa.statistics         #SSR预测分析软件MISA输出的统计文件
|-- SNP_Analysis   #SNP分析结果目录[1]
|   |-- AllSample.SNP_density.png      #SNP密度分布图（png格式）
|   |-- AllSample.SNP_density.stat     #SNP密度分布统计表
|   |-- AllSample.snp.stat             #SNP数量统计表
|   |-- final.indel.vcf                #INDEL提取文件（vcf格式）
|   |-- final.snp.list                 #SNP文件，最主要的结果文件
|   |-- final.snp.vcf                  #SNP提取文件（vcf格式）
|   |-- Pairwised_SNP                  #两两样品间分型不同的SNP目录
|   |   |-- parwised_snp.stat                              #两两样品间分型不同的SNP统计表
|   |   |-- Sample1.Sample2.parwised_snp.hete_hete.list    #样品1和样品2间分型不同的SNP，两样品中均杂合
|   |   |-- Sample1.Sample2.parwised_snp.hete_homo.list    #样品1和样品2间分型不同的SNP，样品1中杂合，样品2中纯合
|   |   |-- Sample1.Sample2.parwised_snp.homo_hete.list    #样品1和样品2间分型不同的SNP，样品1中纯合，样品2中杂合
|   |   |-- Sample1.Sample2.parwised_snp.homo_homo.list    #样品1和样品2间分型不同的SNP，两样品中均纯合
|   |   |-- Sample1.Sample2.parwised_snp.list              #样品1和样品2间分型不同的SNP，所有
|   |   `-- ... ...
|   |-- Sample1.snp.list               #样品1的SNP列表
|   |-- Sample2.snp.list               #样品2的SNP列表
|   `-- ... ...
|-- geneExpression #基因表达量分析结果目录
|   |-- Sample1.Mapped.stat.xls        #样品1测序数据与Unigene或转录本序列的比对统计表
|   |-- Sample1.express.gene_tag.png   #样品1测序数据饱和度模拟图
|   |-- Sample1.fpkm_density.png       #样品1的FPKM密度分布图
|   |-- Sample1.geneExpression.xls     #样品1基因表达量分析结果文件
|   |-- Sample1.insertSize.png         #样品1插入片段长度模拟分布图（svg转png）
|   |-- Sample1.insertSize.r.png       #样品1插入片段长度模拟分布图（ggplot2风格）
|   |-- Sample1.isoforms.results       #样品1转录本表达量分析结果文件
|   |-- Sample1.randcheck.png          #样品1测序Reads在转录本上的位置分布图
|   |-- ... ...
|   |-- Total.gene_tag.png             #各样品测序数据饱和度模拟整合图
|   |-- Total.randcheck.png            #各样品测序Reads在转录本上的位置分布整合图
|   |-- all.fpkm_density.png           #各样品FPKM密度分布对比图
|   |-- all.fpkm_box.png               #各样品FPKM箱线图
|   |-- cor_plot                       #两两样品的表达量相关性散点图目录
|   |   |-- Sample1_vs_Sample2.cor.png           #样品1和样品2的表达量相关性散点图
|   |   `-- ... ...
|   |-- free_com.cor                   #两两样品的表达量相关性（皮尔逊相关系数的平方）统计表
|   `-- sample_cluster.png             #两两样品的表达量相关性热图
|-- DEG_Analysis   #差异表达分析结果目录[2]
|   |-- All_gene_counts.list           #所有基因的表达量（比对片段数）矩阵文件
|   |-- All_gene_fpkm.list             #所有基因的表达量（FPKM值）矩阵文件
|   |-- Condition1_vs_Condition2       #条件1样品和条件2样品的差异表达分析目录[3][4]
|   |   |-- Condition1_vs_Condition2.DEG_cor.png      #差异表达基因相关性散点图
|   |   |-- Condition1_vs_Condition2.FC_FDR.png       #差异表达基因火山图
|   |   |-- Condition1_vs_Condition2.FC_count.png     #差异表达基因MA图
|   |   |-- Condition1_vs_Condition2.DEG_final.xls    #差异表达基因集，即差异表达分析结果文件
|   |   |-- Condition1_vs_Condition2.annotation.xls   #差异表达基因功能注释整合文件
|   |   |-- Condition1_vs_Condition2.cluster.png      #差异表达基因聚类热图
|   |   |-- Condition1_vs_Condition2.cluster.txt      #差异表达基因聚类树形结构（）
|   |   |-- Cog_Anno                                  #差异表达基因COG注释目录
|   |   |   |-- Condition1_vs_Condition2.Cog.classfy.png        #差异表达基因COG注释分类统计图
|   |   |   `-- Condition1_vs_Condition2.Cog.classfy.stat       #差异表达基因COG注释分类统计表
|   |   |-- Kog_Anno                                  #差异表达基因KOG注释目录
|   |   |   |-- Condition1_vs_Condition2.Kog.classfy.png        #差异表达基因KOG注释分类统计图
|   |   |   `-- Condition1_vs_Condition2.Kog.classfy.stat       #差异表达基因KOG注释分类统计表
|   |   |-- Graph                                     #差异表达基因功能富集统计作图目录
|   |   |   |-- Condition1_vs_Condition2.KEGG.Phase.png         #差异表达基因KEGG通路富集散点图
|   |   |   |-- Condition1_vs_Condition2.KEGG.list              #差异表达基因的KEGG富集结果文件
|   |   |   |-- Condition1_vs_Condition2.topGO_BP.pdf           #差异表达基因topGO富集有向无环图（pdf格式，生物学过程）
|   |   |   |-- Condition1_vs_Condition2.topGO_BP.png           #差异表达基因topGO富集有向无环图（png格式，生物学过程）
|   |   |   |-- Condition1_vs_Condition2.topGO_BP.xls           #差异表达基因topGO富集结果文件（生物学过程）
|   |   |   |-- Condition1_vs_Condition2.topGO_CC.pdf           #差异表达基因topGO富集有向无环图（pdf格式，细胞组分）
|   |   |   |-- Condition1_vs_Condition2.topGO_CC.png           #差异表达基因topGO富集有向无环图（png格式，细胞组分）
|   |   |   |-- Condition1_vs_Condition2.topGO_CC.xls           #差异表达基因topGO富集结果文件（细胞组分）
|   |   |   |-- Condition1_vs_Condition2.topGO_MF.pdf           #差异表达基因topGO富集有向无环图（pdf格式，分子功能）
|   |   |   |-- Condition1_vs_Condition2.topGO_MF.png           #差异表达基因topGO富集有向无环图（png格式，分子功能）
|   |   |   `-- Condition1_vs_Condition2.topGO_MF.xls           #差异表达基因topGO富集结果文件（分子功能）
|   |   |-- go_enrichment                             #差异表达分析GO注释富集分析结果目录
|   |   |   |-- Condition1_vs_Condition2.GO.Biological.stat     #差异表达基因Fisher精确检验GO富集结果文件（生物学过程）
|   |   |   |-- Condition1_vs_Condition2.GO.Biological.xls      #差异表达基因Fisher精确检验GO富集统计文件（生物学过程）
|   |   |   |-- Condition1_vs_Condition2.GO.Cellular.stat       #差异表达基因Fisher精确检验GO富集结果文件（细胞组分）
|   |   |   |-- Condition1_vs_Condition2.GO.Cellular.xls        #差异表达基因Fisher精确检验GO富集统计文件（细胞组分）
|   |   |   |-- Condition1_vs_Condition2.GO.Molecular.stat      #差异表达基因Fisher精确检验GO富集结果文件（分子功能）
|   |   |   |-- Condition1_vs_Condition2.GO.Molecular.xls       #差异表达基因Fisher精确检验GO富集统计文件（分子功能）
|   |   |   |-- Condition1_vs_Condition2.GO.list.txt            #差异表达基因编号与注释到的GO节点编号对应表
|   |   |   |-- Condition1_vs_Condition2.GO.png                 #差异表达基因GO二级节点注释统计图（png格式）
|   |   |   |-- Condition1_vs_Condition2.GO.svg                 #差异表达基因GO二级节点注释统计图（svg格式）
|   |   |   `-- Condition1_vs_Condition2.GO_enrichment.stat.xls #差异表达基因GO二级节点注释统计表
|   |   `-- pathway                                   #差异表达基因KEGG功能富集分析结果目录
|   |       |-- kegg_enrichment                                 #差异表达基因KEGG功能富集分析结果目录
|   |       |   |-- Condition1_vs_Condition2.KEGG.png           #差异表达基因KEGG分类图（png格式）
|   |       |   |-- Condition1_vs_Condition2.KEGG.stat          #差异表达基因KEGG分类统计表
|   |       |   |-- Condition1_vs_Condition2.KEGG.svg           #差异表达基因KEGG分类图（svg格式）
|   |       |   |-- Condition1_vs_Condition2.KEGG.tree.stat     #差异表达基因KEGG分类统计详表
|   |       |   |-- Condition1_vs_Condition2.KEGG.xls           #KEGG通路注释结果统计文件
|   |       |   `-- Condition1_vs_Condition2.Kegg.ko            #差异表达基因KEGG注释结果文件
|   |       `-- kegg_map                                        #差异表达基因注释到的KEGG通路图目录
|   |           |-- ko00010.html                                     #差异表达基因注释到的KEGG通路图（html格式）
|   |           |-- ko00010.png                                      #差异表达基因注释到的KEGG通路图（png格式）
|   |           `-- ... ...
|   |-- ... ...
|   |-- All_DEG                        #各差异表达基因集并集分析目录
|   |   |-- All.DEG_final.xls          #各差异表达基因集并集基因的表达量文件
|   |   |-- All.DEG.cluster.png        #各差异表达基因集并集基因的表达聚类热图      
|   |   |-- All.DEG.cluster.txt        #各差异表达基因集并集基因的表达聚类树形图       
|   |-- All_DEG_veen.genes             #各差异表达基因集元素列表[5]
|   |-- All_DEG_veen.png               #各差异表达基因集韦恩图[5]
|   |-- All_DEG_veen.stat              #各差异表达基因集数目统计表[5]
|   |-- DEG.anno.stat                  #各差异表达基因集的功能注释统计表
|   `-- DEG.stat                       #各差异表达基因集数目统计表
`-- readme.txt     #说明文档

注：
[1] 针对无参考转录组分析，只有样本个体数目在两个及以上时才能进行SNP分析，因此，只有项目有两个以上的样品时才有SNP分析结果，否则不会生成该目录；
[2] 同样，差异表达分析的是两个条件下（即不同时间、组织或处理等）样本间基因表达水平的差异及显著性，因此必须要求至少有两个样品，否则不会生成该目录；
[3] 对于设立生物学重复（Biological Replicates）的实验，同一条件会有多个样品（通常3个及以上），我们采用DESeq软件进行两条件样品组间的差异表达分析；而对于没有设立生物学重复的实验，同一条件只有1个样品，我们采用EBSeq软件进行两个样品间的差异表达分析。
[4] 在分析结果中，使用“A_vs_B”的方式命名差异表达基因集，如T1_vs_T2或T1_T2_vs_T3_T4等。通常情况下，对于两个样品之间的差异表达基因集，A表示对照样品、野生型样品或前一个时间点样品；而B表示对应的处理样品、突变型样品或后一个时间点样品。相应地，对于两个条件（即两组样品）之间的差异表达基因集，A表达含有多个重复样品（Duplicates）的对照组、野生型组或前一个时间点样品组；B表示对应的处理组、突变型组、后一个时间点样品组。根据两（组）样品之间表达水平的相对高低，差异表达基因可以划分为上调基因（Up-regulated Gene）和下调基因（Down-regulated Gene）。样品（组）B中的表达水平高于样品（组）A中的基因称之为上调基因；反之为下调基因。因此，上调和下调是相对的，由所给A和B的顺序决定，若更换A和B的顺序会完全反过来，但这不会对分析结果产生实质性的影响。
[5] 仅当差异表达基因集数目在2到5个之间才能绘制韦恩图，才生产该文件。


基本文件格式说明：
********************************************************************************
1. FASTQ格式：
=================================================
FASTQ格式是一种用于存放核酸序列及其对应的质量值的文本格式，是存储高通量测序数据的标准格式。我们的测序数据均以FASTQ格式存储，通常以“.fastq” 或 “.fq”为文件后缀（以“.fq.gz”为文件后缀的是压缩过的FASTQ格式文件）。格式如下：

@HWI-7001455:116:H8E71ADXX:2:2212:8408:20433 1:Y:0:ATCACG
CTTCAACCAGGTCACCGGCATCAACGTCATCAACTTCTACGCGCCGTTCATGTTCCGGACCATCGGGCTCAAGGAGAGCGCGTCCCTCATGTCGGCCGTGG
+
;<5;=<??#############################################################################################
@HWI-7001455:116:H8E71ADXX:2:2107:14756:56485 1:Y:0:ATCACG
CAGGTGCTGACAGCAATCGGTAACTTCAGCATCTGCTCCATTGGGGTCGGCATCCTCGTCGAGATCATCGTCATGTTCCCAATCCAGCACCGGAAGTACCG
+
;5;?#################################################################################################

FASTQ文件中通常每4行对应一个序列单元：
第一行 以@开头，后面接着序列标识（ID）以及其他可选的描述信息；
第二行 序列行；
第三行 以“+”开头，后面接着可选的描述信息；
第四行 第二行每个字母对应的质量值编码，长度必须和第二行的序列长度相同。


2. FASTA格式：
=================================================
FASTA格式是一种用于表示核苷酸序列或多肽序列的标准文本格式，其中核苷酸或氨基酸用单个字母表示，序列之前放置序列名字和注释信息。我们组装得到的Contig序列和转录本序列、构建的Unigene序列、编码区序列（CDS）以及对应的氨基酸序列等都以FASTA格式存储。FASTA格式文件通常以“.fasta” 或 “.fa”结尾。以CDS文件为例：

>c1000.graph_c0|orf1 type=3prime_partial len=462nt loc=c1000.graph_c0:62-523:-
TTGGACAGGTCTTGGATGGTTTTGGAGACTGCAGGGGGTCACGGCCGCAACCCCAATGGC
ATCCTCGGCCTTCTATGCCGGGAACACTTCCCTGGGCTTGTCGAGTACGCCGGAGTGACG
AGCCCAGCGTACACCTTCGACCACTACGCCGTCGCCCCCGATGCAGTAGATCGGGACGGC
AGACAATTCAACAACAAGGCGGAGCGGGTCAAGCAAGAGCTGTGGGTAAGTCTTCCTCGC
ACTATATTGTTTAATAAGTCGCATTTGTTGCATATTCTTGAAATAATGTATGGATACATC
GTCTTTGTATGTAGGATTTCTTCAGGTGCGATGCTGGATACGAGGCCAGGGCGGATGTGG
TGTCTACGACGTGCTGTAAGAAGCTCGTCGTGGACATGCACTACGAGGCGCGCATCCAGG
CCATCGTCACTTACCACGGCTCCGTCCTTGGGGAGAAGGTGA
>c1020.graph_c0|orf1 type=complete len=705nt loc=c1020.graph_c0:550-1254:+
ATGCTGAATCTAAGTAGCTCTGCTGTATCTGCGCCATCAAGAACTCATGTGGATCATGCC
GCTCTTACCGGTGCCTCTCATCCAGCTTCTACGGTCAAGACACGCATGTGCACCAAGTAC
AACACTACAGAAGGCTGCAAGTTCGGTGATAAGTGCCATTTCGCTCATAGCGAAAGAGAG
CTTGCGAAGCCAGCCTACATGTCTCAAGAAGGACCTCCTATGGGTGGTCGATATGGACGA
GCTGAACCTATGCAACCAGCTGCCATGGGCCCTCCAGCAGGAAACTTTGGTGCCTCGGCG
ACTGCCAAGATCAGTGTGGACGCCTCTCTGGCCGGTGGCATAATCGGCAAGGGTGGGGTC
AACACTAAGCAGATAAGTAGAATTACAGGCGTCAAGCTCTCCATCCGCGACCACGAGTCT
AACCCCAACCTAAAGAACATCGAGCTGGAAGGCAATTTTGACCAGATCAAGCAAGCCAGC
GACTTGGTGCGTGATCTCATCGCAACCATCAGCGCAAGCATGCCAGCGAAAAATCCATCT
GCTGCCGCGGCACCAGCAGGAGGAGGCCGAGGTGGTGGTCCAGGGGGCAAGAGCAGCAAC
TACAAGACGAAGCTGTGCGAGAACTTCTTGAAGGGTGCCTGTACTTTTGGTGACCGGTGC
CACTTCGCCCATGGCGAGACGGAGCAGCGGAGAGGTGCTGCATGA

FASTA格式每一个序列单元由描述信息和序列数据两部分构成，每一个序列单元从“>”（大于号）开始到下一个“>”之前结束。
描述信息与“>”放于第一行，“>”后面紧接序列标识（ID，如c1000.graph_c0|orf1），二者之间不能有空格，之后是可选的描述；
序列数据从第二行开始，可放于单行或多行，通常每一行为60或100个字符，直到出现“>”之前为止。
“>”标志着一个新的序列单元的开始。


重要文本文件说明：
********************************************************************************
1. 测序数据评估统计表【/cleandata/AllSample_GC_Q.stat】各列说明：
=================================================
SampleID：样品信息单样品名称；
ReadSum：Clean Data中pair-end Reads总数；
BaseSum：Clean Data总碱基数；
GC(%)：Clean Data GC含量，即Clean Data中G和C两种碱基占总碱基的百分比；
N(%)：Clean Data中未识别碱基占总碱基的百分比；
Q20(%)：Clean Data质量值大于或等于20的碱基所占的百分比；
Q30(%)：Clean Data质量值大于或等于30的碱基所占的百分比。

2. 测序组装结果统计表【/Assembly/All_Combination/*/All_Combination.*.stat.xls】各行说明：
=================================================
Length Range：Contig/Transcript/Unigene序列的不同长度区间；
Total Number：组装得到的Contig/Transcript/Unigene序列的总数；
Total Length：组装得到的Contig/Transcript/Unigene序列的总长度；
N50 Length：Contig/Transcript/Unigene序列的N50的长度；
Mean Length：Contig/Transcript/Unigene序列的平均长度。

3. 基因可变剪接统计表【/Assembly/All_Combination/Unigenes/All_Combination.Unigene.Cluster.stat.xls】各列说明：
=================================================
AS_number：基因检测到的可变剪接（Alternative Splicing Event）数目，即基因所拥有的转录本数目减1；
Gene_number：检测到特定可变剪接数目的基因数目；
Percent：检测到特定可变剪接数目的基因在总基因数目中所占的百分比。

4. 功能注释整合文件【/Unigene/Unigene_Anno/Integrated_Function.annotation.xls】各列（可选）说明：
=================================================
GeneID：获得功能注释的Unigene编号；
COG_class：COG库注释分类；
COG_class_annotation：COG库注释信息；
KOG_class：KOG库注释分类；
KOG_class_annotation：KOG库注释信息；
GO_annotation：GO数据库注释信息；
KEGG_annotation：KEGG库注释信息；
Pfam_annotation：Pfam库注释信息；
Swissprot_annotation：Swiss-Prot库注释信息；
nr_annotation：nr库注释信息；
TrEMBL_annotation：TrEMBL库注释信息；
nt_annotation：nt库注释信息。

5. 注释数目统计表【/Unigene/Unigene_Anno/Function_Annotation.stat.xls】各列说明：
=================================================
Anno_Database：用于注释的数据库；
Annotated_Number：获得对应数据库注释信息的Unigene总数目；
300<=length<1000：获得对应数据库注释信息且长度在300nt和1000nt之间的Unigene数目；
length>=1000：获得对应数据库注释信息且长度1000nt即以上的Unigene数目。

6. COG或KOG数据库注释结果文件【/Unigene/Unigene_Anno/Maize.Unigene.fa.[C|K]og_class.txt】各列说明：
=================================================
Gene name：获得COG或KOG数据库注释信息的Unigene编号；
Portein name in COG/KOG：COG或KOG数据库中Unigene比对到的蛋白名称；
E_value：BLAST软件比对E值，表示其它序列与目标序列（query，即Unigene序列）相似度要大于这条显示的序列（subject，即数据库中Unigene比对到的蛋白序列）的可能性，值越低越好；
Identity：BLAST软件比对Identity，即比对上的碱基中完全匹配的碱基比例；
Score：BLAST软件比对打分；
Organism：COG或KOG数据库中Unigene比对到的蛋白物种来源；
COG/KOG id：COG或KOG数据库中Unigene比对到的蛋白编号；
COG/KOG class defination：COG或KOG分类定义；
Function code：分类缩写；
Functional categories：功能类别；
Function class defination：功能分类定义。

7. COG或KOG注释分类统计表【/Unigene/Unigene_Anno/SpeciesName.Unigene.fa.[C|K]og.cluster.stat】对应COG或KOG注释分类统计图。各列说明：
=================================================
ID：COG或KOG注释分类名称缩写；
Class_Name：COG或KOG注释分类名称；
Numbers：注释到对应分类上的Unigene数目。

8. GO数据库注释结果文件【/Unigene/Unigene_Anno/SpeciesName.Unigene.fa.GO.anno.txt】各列说明：
Gene：Unigene编号；
GO_Anno：GO注释信息，即Unigene注释到的GO节点（包括节点分类、语义描述和节点编号），之间用“;”分隔。

9. GO二级节点注释统计表【/Unigene/Unigene_Anno/SpeciesName.Unigene.fa.GO_enrichment.stat.xls】对应GO二级节点注释统计图，各列说明：
=================================================
GO_classify：GO二级节点语义或名称；
SpeciesName.Unigene：注释到给二级节点的Unigene数目。

10.KEGG数据库注释结果文件【/Unigene/Unigene_Anno/SpeciesName.Unigene.fa.Kegg.ko】各列说明：
=================================================
Gene_id：Unigene编号；
第二列用“|”分隔的依次为：
KO：KEGG数据库中的序列K系统编号；
e_value：BLAST软件比对E值，表示其它序列与目标序列（query，即Unigene序列）相似度要大于这条显示的序列（subject，即数据库中Unigene比对到的蛋白序列）的可能性，值越低越好；
Database_Genes：KEGG数据库中的序列编号，用“:”连接物种名称缩写和基因编号而成；
Anno：KEGG数据库中的序列注释信息。

11. KEGG通路与注释到该通路的Unigene编号对应表【/Unigene/Unigene_Anno/SpeciesName.Unigene.fa.Kegg.pathway】各列说明：
=================================================
pathway：KEGG通路名称；
pathway_id：KEGG通路编号；
Gene_number：注释到对应通路上的Unigene编号，之间用“;”分隔；
Gene_id KOs：Unigene对应的KEGG库序列的K系统编号，之间用“+”分隔。

12. Pfam数据库注释结果文件【/Unigene/Unigene_Anno/SpeciesName.Unigene.fa.Pfam.anno.txt】各列说明：
=================================================
Gene_ID：Unigene编号；
Pfam_IDs：Pfam家族编号；
Pfam_Description：Pfam家族名称。

13. Swiss-Prot数据库注释结果文件【/Unigene/Unigene_Anno/SpeciesName.Unigene.fa.Swissprot.anno.txt】各列说明：
=================================================
SwissprotGeneID：Unigene编号；
Database_ID：KEGG数据库中Unigene比对到的蛋白编号；
E_value：BLAST软件比对E值，表示其它序列与目标序列（query，即Unigene序列）相似度要大于这条显示的序列（subject，即数据库中Unigene比对到的蛋白序列）的可能性，值越低越好；
Identity：BLAST软件比对Identity，即比对上的碱基中完全匹配的碱基比例；
Score：BLAST软件比对打分；
Annotation：Swiss-Prot数据库中的序列注释信息。

14. nr数据库注释结果文件【/Unigene/Unigene_Anno/SpeciesName.Unigene.fa.nr.anno.txt】各列说明：
=================================================
NrGeneID：Unigene编号；
Database_ID：KEGG数据库中Unigene比对到的蛋白编号；
E_value：BLAST软件比对E值，表示其它序列与目标序列（query，即Unigene序列）相似度要大于这条显示的序列（subject，即数据库中Unigene比对到的蛋白序列）的可能性，值越低越好；
Identity：BLAST软件比对Identity，即比对上的碱基中完全匹配的碱基比例；
Score：BLAST软件比对打分；
Annotation：nr数据库中的序列注释信息。

15. SSR分析结果整合表【/Unigene/Unigene_SSR/SpeciesName.Unigene.1000.fa.SSR.result.xls】各列说明：
=================================================
Gene_ID：Unigene编号；
SSR_nr：同一Unigene上的SSR序号；
SSR_type：SSR类型，包括完美单碱基重复（p1）、完美双碱基重复（p2）、完美三碱基重复（p3）、完美四碱基重复（p4）、完美五碱基重复（p5）、完美六碱基重复（p6）和混合SSR（c，即包含至少两个完美SSR，且之间距离小于100bp）；
SSR：SSR序列，括号内为重复单元，括号外数字表示重复次数；
SSR_Start：SSR在Unigene上的开始位置；
SSR_End：SSR在Unigene上的结束位置。
FPr1(5'-3')：第一条正向引物序列；
Tm：第一条正向引物序列的退火温度，单位为°C；
Size：第一条正向引物序列的长度；
RPr1(5'-3')：第一条反向引物序列；
Tm：第一条反向引物序列的退火温度，单位为°C；
Size：第一条反向引物序列的长度；
PSize：产物的长度；
Start：产物在基因上的开始位置；
End：产物在基因上的结束位置。
FPr2(5'-3')：第二条正向引物序列；
... ...

16. 样品1的SNP列表【/SNP_Analysis/final.snp.list】各列说明：
=================================================
GeneID：Unigene编号；
Pos：SNP位点在Unigene上的位置；
Ref：Unigene上的SNP等位；
Alt：测序样品中识别到的其他的SNP等位；
T1：样品1该SNP位点的分型；
Depth：样品1该SNP位点的测序深度；
AlleDepth：样品1该SNP位点的各等位测序深度，用“,”分隔。
T2：样品2该SNP位点的分型；
Depth：样品2该SNP位点的测序深度；
AlleDepth：样品2该SNP位点的各等位测序深度，用“,”分隔。
... ...
Effect：（缺失信息）
Codon_change：（缺失信息）
Gene_id：（缺失信息）

17. 样品1和样品2间分型不同的SNP【/SNP_Analysis/Pairwised_SNP/Sample1.Sample2.parwised_snp.list】各列说明：
=================================================
GeneID：Unigene编号；
Postion：SNP位点在Unigene上的位置；
RefAlle：Unigene上的SNP等位；
AltAlles：测序样品中识别到的其他的SNP等位；
Sample1.genotype：样品1该SNP位点的分型；
Sample1.totalDep：样品1该SNP位点的测序深度；
Sample1.alleDeps：样品1该SNP位点的各等位测序深度，用“,”分隔。
Sample2.genotype：样品2该SNP位点的分型；
Sample2.totalDep：样品2该SNP位点的测序深度；
Sample2.alleDeps：样品2该SNP位点的各等位测序深度，用“,”分隔。

18. 样品1基因表达量分析结果文件【/geneExpression/Sample1.geneExpression.xls】各列说明：
=================================================
Gene_ID：Unigene编号；
Length：Unigene的长度；
Effective_Length：Unigene有效长度，即该基因不同转录本的平均长度；
TPM：TPM方法标准化后的基因表达丰度值；
FPKM：FPKM方法标准化后的基因表达丰度值；
Transcript_ID(s)：转录本的编号；
Expected_Count：标准化后的片段数。
======================或=========================
Gene_ID：Unigene的编号；
Length：Unigene的长度；
Depth：测序Reads在Unigene上的平均深度；
Coverage：测序Reads在Unigene上的覆盖度；
FPKM：FPKM方法标准化后的基因表达丰度值；
TotalReads：比对到Unigene上的测序Reads数目；
UniqReads：比对到Unigene唯一位置上的测序Reads数目；
MultiReads：比对到多个Unigene或一个Unigene多个位置上的测序Reads数目。

19. 样品1测序数据与Unigene或转录本序列的比对统计表【/geneExpression/Sample1.Mapped.stat.xls】各行说明：
=================================================
Total Reads：测序Reads总数目；
Mapped Reads：比对到Unigene上的测序Reads总数目，及其在测序Reads总数目中所占的百分比；
Uniq mapped Reads：比对到Unigene唯一位置上的测序Reads总数目，及其在测序Reads总数目中所占的百分比；
Multi mapped Reads：比对到多个Unigene或一个Unigene多个位置上的测序Reads总数目，及其在测序Reads总数目中所占的百分比。

20. 所有基因的表达量（比对片段数）矩阵文件【/DEG_Analysis/All_gene_counts.list】各列说明：
=================================================
ID：Unigene编号；
Sample1：样品1中基因的表达量（比对片段数）；
... ...：其他样品中基因的表达量（比对片段数）；
geneLength：Unigene序列长度。

21. 所有基因的表达量（FPKM值）矩阵文件【/DEG_Analysis/All_gene_fpkm.list】各列说明：
=================================================
ID：Unigene编号；
Sample1：样品1中基因的表达量（FPKM值）；
... ...：其他样品中基因的表达量（FPKM值）。

22. 差异表达基因集，即差异表达分析结果文件【/DEG_Analysis/Condition1_vs_Condition2/Condition1_vs_Condition2.DEG_final.xls】各列说明：
=================================================
GeneID：基因编号；
... ...：对应样品中基因的表达量FPKM值；
FDR：错误发现率；
log2FC：表达量差异倍数的对数值；
regulated：上调基因（up）还是下调基因（down）。

23. 差异表达基因topGO富集结果文件【/DEG_Analysis/Condition1_vs_Condition2/Graph/Condition1_vs_Condition2.topGO_*.xls】各列说明：
=================================================
GO.ID：GO 节点的编号；
Term：GO节点名称；
Annotated：所有基因注释到该功能的基因数；
Significant：DEG注释到该功能的基因数；
Expected：注释到该功能DEG数目的期望值；
KS：富集节点的显著性统计，KS值越小，表明富集越显著。

24. 差异表达基因Fisher精确检验GO富集统计文件【/DEG_Analysis/Condition1_vs_Condition2/go_enrichment/Condition1_vs_Condition2.GO.*.xls】各列说明：
=================================================
Term：GO节点描述和编号；
mNum：注释到该GO节点的差异表达基因数目；
nNum：注释到该GO节点的Unigene数目；
MNum：注释到GO数据库的差异表达基因总数；
NNum：注释到GO数据库的Unigene总数；
GIDs：注释到该GO节点的差异表达基因编号，之间用“;”分隔。

25. 差异表达基因Fisher精确检验GO富集结果文件【/DEG_Analysis/Condition1_vs_Condition2/go_enrichment/Condition1_vs_Condition2.GO.*.stat】各列说明：
=================================================
Gene_Ontology_term：GO节点描述和编号；
Cluter_frequency：注释到该GO节点的差异表达基因数在注释到GO数据库的差异表达基因总数中所占的比例；
Genome_frequency：注释到该GO节点的Unigene数在注释到GO数据库的Unigene总数中所占的比例；
P-value：Fisher精确检验富集显著性P值；
Corrected_P-value：富集显著性P值矫正值。




-------------------------------------------------------------------------------
tree v1.5.0 1996-2009 by Steve Baker
copyright (c) BMK 2014 by Simon Young
last update 2015-04-08
