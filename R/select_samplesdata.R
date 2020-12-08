
getwd()
setwd('C:\Users\leekyle\Downloads\Desktop\sangho')

####################################################conect mysql in r##########################################################

install.packages('RMySQL')
library(RMySQL)
drv <- dbDriver("MySQL")
con <- dbConnect(MySQL(), 
                 dbname = "tcga", 
                 user = "hwanghou", 
                 host = "103.22.220.149", 
                 port = 13306,
                 password = "yonsei2020!")

### tumor samples
r_sql <- c("SELECT `icgc_donor_id`,`project_code`,`submitted_sample_id`,`gene_stable_id`,`normalized_expression_level` FROM `gene_expression` WHERE `submitted_sample_id`like'%-01A%'")
data <- dbGetQuery(con, r_sql) #running 20min

#####################################################################################
#########################          data check           #############################
library(data.table)

head(table(data$submitted_sample_id), 10)
a<-as.data.table(table(data$submitted_sample_id))
summary(a$N)
table(a$N)
nrow(a)



# to much gene (>2000)
gene_2000up_id<-subset(a, N>=20000) 
nrow(gene_2000up_id) #832/8775
sample<-subset(data, submitted_sample_id=="TCGA-05-4244-01A-01R-1107-07")
nrow(sample) #34959
head(sample)
sample<-sample[c(order(sample$gene_stable_id)),]
head(sample, 100)
summary((as.data.table(table(sample$gene_stable_id))$N))
subset((as.data.table(table(sample$gene_stable_id))), N==4) #1
subset((as.data.table(table(sample$gene_stable_id))), N==3) #5
str(subset((as.data.table(table(sample$gene_stable_id))), N==2)) #49
summary((as.data.table(table(sample$gene_stable_id))$N))
str(unique(sample$gene_stable_id)) #34897

# to much gene (<20000)
gene_2000down_id<-subset(a, N<=20000) #832/8775
nrow(gene_2000down_id) #7944/8775
head(gene_2000down_id)
#2
sample2<-subset(data, submitted_sample_id=="TCGA-02-0009-01A-01R-0177-01")
nrow(sample2) #11968
str(unique(sample2$gene_stable_id)) #11922
a<-rbind(as.data.table(unique(sample2$gene_stable_id)), as.data.table(unique(sample$gene_stable_id)))
nrow(a) # 34897+11922
a2<-as.data.table(table(a$V1))
head(a2, 50)
summary(a2$N)
a3<-subset(a2, N==2)
a4<-subset(a2, N==1)
write.csv(a3, file = "a3.csv", row.names = FALSE)
write.csv(a4, file = "a4.csv", row.names = FALSE)

#3 (no gene %ENSG%)
tail(gene_2000down_id, 20)
sample3<-subset(data, submitted_sample_id=="TCGA-P5-A5F1-01A-11R-A28M-07")
head(sample3)
nrow(sample3) #18457
str(unique(sample3$gene_stable_id)) #18456
sample3[grep("ENSG", sample3$gene_stable_id),]

#####################################################################################
######################### match gene id with gene name  #############################

str(unique(data$submitted_sample_id))
str(unique(stage$bcr_patient_barcode))

##### exclude "ENSG"
library(dplyr)
sample_ENSG_na<-data %>% 
  filter(!grepl('ENSG', gene_stable_id))
head(sample_ENSG_na, 20)
sample_ENSG_na[grep("ENSG", sample_ENSG_na$gene_stable_id),]

str(unique(sample_ENSG_na$submitted_sample_id)) #5575/8775
str(unique(sample_ENSG_na$gene_stable_id)) #20458

x<-as.data.table(table(sample_ENSG_na$submitted_sample_id))
summary(x$N)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#18457   18457   18457   18570   18457   36914
table(x$N)
#18457 21131 36914
#5492    57    26

y<-as.data.table(table(sample_ENSG_na$gene_stable_id))
head(y)
summary(y$N)
table(y$N) #5575+26 = 5601 and 18294+'SLC35E2'= 18295
#  57  5544  5601 11316('SLC35E2') 47424('?') 
#2001   161 18294     1                1 

sample_ENSG_na_genelist<-subset(y, N==5601 | N==11316)
nrow(sample_ENSG_na_genelist) #gene list 18295

##############################
###check gene 18457 sampels
x18457<-subset(x, N==18457)

head(x18457)
colnames(x18457) = c("submitted_sample_id", "gene_n")
head(x18457)
x18457_allsamples<-merge(x=x18457, y=data, by = "submitted_sample_id", all.x = TRUE)
nrow(x18457_allsamples) # 113939354
x18457_allsamples2<-x18457_allsamples %>% 
  filter(!grepl('ENSG', gene_stable_id))
x18457_allsamples2[grep("ENSG", x18457_allsamples2$gene_stable_id),]

nrow(x18457_allsamples2) # 18457(gene)*4838(samples) = 101365844
str(unique(x18457_allsamples2$gene_stable_id)) #18456 ("SLC35E2" expression data is 2)
d<-as.data.table(table(x18457_allsamples2$gene_stable_id))
head(d)
table(d$N)
subset(d, N==10984) 
#        V1     N
#1: SLC35E2 10984


##############################
### check gene 21131 sampels
x21131<-subset(x, N==21131)
nrow(x21131) #samples = 57
head(x21131)
colnames(x21131) = c("submitted_sample_id", "gene_n")
head(x21131)

x21131_allsamples<-merge(x=x21131, y=data, by = "submitted_sample_id", all.x = TRUE)
nrow(x21131_allsamples) # 1204467
x21131_allsamples2<-x21131_allsamples %>% 
  filter(!grepl('ENSG', gene_stable_id))
x21131_allsamples2[grep("ENSG", x21131_allsamples2$gene_stable_id),]

nrow(x21131_allsamples2) # 21131(gene)*57(samples) = 1204467
str(unique(x21131_allsamples2$gene_stable_id)) # 20297 gene
d<-as.data.table(table(x21131_allsamples2$gene_stable_id))
head(d)
table(d$N)
#   57   228 47424 
#20295     1     1

subset(d, N==228) # ("SLC35E2" expression data is 4)
#        V1   N
#1: SLC35E2 228
subset(d, N==47424) # ("?" expression data is 4)
#   V1     N
#1:  ? 47424

##############################
### check gene 36914 sampels
x36914<-subset(x, N==36914)
nrow(x36914) #samples = 26
head(x36914)
colnames(x36914) = c("submitted_sample_id", "gene_n")
head(x36914)
x36914_allsamples<-merge(x=x36914, y=data, by = "submitted_sample_id", all.x = TRUE)
nrow(x36914_allsamples) # 959764
x36914_allsamples2<-x36914_allsamples %>% 
  filter(!grepl('ENSG', gene_stable_id))
x36914_allsamples2[grep("ENSG", x36914_allsamples2$gene_stable_id),]

nrow(x36914_allsamples2) # 36914(gene)*26(samples) = 959764
str(unique(x36914_allsamples2$gene_stable_id)) # 18456 gene
d<-as.data.table(table(x36914_allsamples2$gene_stable_id))
head(d)
table(d$N)
#   52   104 
#18455     1 

subset(d, N==104) # ("SLC35E2" expression data is 2)
#        V1   N
#1: SLC35E2 104


#####   #####   #####  #####  
#####  include "ENSG"  ##### 
#####   #####   #####  #####  
sample_ENSG<-data[grep("ENSG", data$gene_stable_id),]
sample_ENSG %>% 
  filter(!grepl('ENSG', gene_stable_id))
str(unique(sample_ENSG$submitted_sample_id)) #3948/8775
str(unique(sample_ENSG$gene_stable_id)) #17432

head(sample_ENSG)

x<-as.data.table(table(sample_ENSG$submitted_sample_id))
summary(x$N)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#11968   11968   16813   15563   17321   17321
table(x$N)
#11968 16502 16595 16702 16714 16717 16718 16730 16742 16760 16762 16764 16766 16767 16771 16772 16775 16776 16778 16780 16781 16782 16783 16786 16787 16788 16789 16790 16791 16792 
# 1124     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     2     2     4     1     3     3     1     1     2     2     5     2 
#16793 16794 16795 16796 16797 16798 16799 16800 16801 16802 16803 16804 16805 16806 16807 16808 16809 16810 16811 16812 16813 17321 
#    2     1     3     7     4     6     7     1    10     6    11    10    18    27    23    40    50    87   134   276  1039  1017
str(unique(sample_ENSG$gene_stable_id)) #17432 gene

y<-as.data.table(table(sample_ENSG$gene_stable_id))
head(y)
summary(y$N)
table(y$N)
#1017  1124  1770  1789  1797  1800  1803  1804  1805  1806  1807  2034  2141  2248  2738  2758  2782  2783  2786  2788  2789  2792  2794  2795  2796  2797 
# 542     5     1     1     2     1     3     2     5     5    94    16   114     1     1     1     1     1     1     1     2     1     1     2     3     6 
#2799  2800  2801  2802  2803  2804  2805  2806  2807  2808  2809  2810  2811  2812  2813  2814  2815  2816  2817  2818  2819  2820  2821  2822  2823  2824 
#   2     3     5     1     1     4     8     4     3     5     3     7     5    10     9     9    10    13     9     8    23    35    37    70   175  4330 
#2929  2930  2931  3051  3614  3889  3900  3907  3910  3915  3917  3920  3921  3924  3925  3927  3928  3929  3930  3931  3932  3933  3934  3935  3936  3937 
#   1     2    66     2     2     1     1     1     1     1     3     1     1     2     1     1     3     1     3     4     3     4     3     5     2     5 
#3938  3939  3940  3941  3942  3943  3944  3945  3946  3947  3948  5421  5602  5644  5646  5648  5860  7896  8793 11844 15792 25692 
#   7     2     3     5     5    23    18    36    52   144 11361     2     1     1     1    19     1    24     1     2     1     1 



##############################
### check gene 11968 sampels
x11968<-subset(x, N==11968)
nrow(x11968) #samples = 1124
head(x11968)
colnames(x11968) = c("submitted_sample_id", "gene_n")
head(x11968)

x11968_allsamples<-merge(x=x11968, y=data, by = "submitted_sample_id", all.x = TRUE)
nrow(x11968_allsamples) # 13452032
x11968_allsamples2<-x11968_allsamples[grep("ENSG", x11968_allsamples$gene_stable_id),]
x11968_allsamples2 %>% 
  filter(!grepl('ENSG', gene_stable_id))
nrow(x11968_allsamples2) # 11968(gene)*1124(samples) = 13452032
str(unique(x11968_allsamples2$gene_stable_id)) # 11922 gene
d<-as.data.table(table(x11968_allsamples2$gene_stable_id))
nrow(d)
head(d, 100)
table(d$N)
# 1124  2248  3372  4496 13488 
#11891    26     3     1     1
subset(d, N>=2248)
x11968_dupli<-subset(d, N>=2248) 

#               V1     N
#1: ENSG00000015568  3372
#2: ENSG00000059145  2248
#3: ENSG00000072310  2248
#4: ENSG00000154864  4496
#5: ENSG00000168477  2248
#6: ENSG00000172352  2248
#7: ENSG00000184779  2248
#8: ENSG00000185565  2248
#9: ENSG00000187775  3372
#10: ENSG00000189013  2248
#11: ENSG00000197932  3372
#12: ENSG00000198040  2248
#13: ENSG00000198393  2248
#14: ENSG00000204344  2248
#15: ENSG00000204501  2248
#16: ENSG00000206258  2248
#17: ENSG00000206342  2248
#18: ENSG00000215790  2248
#19: ENSG00000226033  2248
#20: ENSG00000226257  2248
#21: ENSG00000229341  2248
#22: ENSG00000229353  2248
#23: ENSG00000231608  2248
#24: ENSG00000233323  2248
#25: ENSG00000234947  2248
#26: ENSG00000235173  2248
#27: ENSG00000236221  2248
#28: ENSG00000236236  2248
#29: ENSG00000236250  2248
#30: ENSG00000241476  2248
#31: ENSG00000250589 13488


##############################
### check gene 16813 sampels

x16813<-subset(x, N==16813)
nrow(x16813) #samples = 1039
head(x16813)
colnames(x16813) = c("submitted_sample_id", "gene_n")
head(x16813)
x16813_allsamples<-merge(x=x16813, y=data, by = "submitted_sample_id", all.x = TRUE)
nrow(x16813_allsamples) # 24039399
x16813_allsamples2<-x16813_allsamples[grep("ENSG", x16813_allsamples$gene_stable_id),]
x16813_allsamples2 %>% 
  filter(!grepl('ENSG', gene_stable_id))
nrow(x16813_allsamples2) # 16813(gene)*1039(samples) = 17468707
str(unique(x16813_allsamples2$gene_stable_id)) # 16751 gene
d<-as.data.table(table(x16813_allsamples2$gene_stable_id))
nrow(d)
head(d, 100)
table(d$N)
# 1039  2078  3117  4156 
#16696    49     5     1
subset(d, N>=2078) 
x16813_dupli<-subset(d, N>=2078) 

'''
             V1    N
1: ENSG00000015568 3117
2: ENSG00000059145 2078
3: ENSG00000072310 2078
4: ENSG00000109181 2078
5: ENSG00000154864 4156
6: ENSG00000164871 2078
7: ENSG00000168477 2078
8: ENSG00000171489 2078
9: ENSG00000172058 2078
10: ENSG00000172352 2078
11: ENSG00000182330 3117
12: ENSG00000183385 2078
13: ENSG00000183558 2078
14: ENSG00000184779 2078
15: ENSG00000185565 2078
16: ENSG00000185684 2078
17: ENSG00000187243 2078
18: ENSG00000187627 2078
19: ENSG00000187653 2078
20: ENSG00000187775 3117
21: ENSG00000189013 2078
22: ENSG00000197932 3117
23: ENSG00000198040 2078
24: ENSG00000198393 2078
25: ENSG00000198457 2078
26: ENSG00000203811 2078
27: ENSG00000203989 2078
28: ENSG00000204164 2078
29: ENSG00000204344 2078
30: ENSG00000204828 2078
31: ENSG00000205176 2078
32: ENSG00000205944 2078
33: ENSG00000206258 2078
34: ENSG00000206338 2078
35: ENSG00000206342 2078
36: ENSG00000215790 2078
37: ENSG00000225932 2078
38: ENSG00000226033 2078
39: ENSG00000226257 2078
40: ENSG00000229341 2078
41: ENSG00000229353 2078
42: ENSG00000231608 2078
43: ENSG00000231852 2078
44: ENSG00000232414 2078
45: ENSG00000233151 2078
46: ENSG00000233323 2078
47: ENSG00000233630 3117
48: ENSG00000234947 2078
49: ENSG00000235134 2078
50: ENSG00000235173 2078
51: ENSG00000236221 2078
52: ENSG00000236236 2078
53: ENSG00000236250 2078
54: ENSG00000240428 2078
55: ENSG00000241476 2078
'''

##############################
### check gene 17321 sampels
x17321<-subset(x, N==17321)
nrow(x17321) #samples = 1017
head(x17321)
colnames(x17321) = c("submitted_sample_id", "gene_n")
head(x17321)
x17321_allsamples<-merge(x=x17321, y=data, by = "submitted_sample_id", all.x = TRUE)
nrow(x17321_allsamples) # 17615457
x17321_allsamples2<-x17321_allsamples[grep("ENSG", x17321_allsamples$gene_stable_id),]
x17321_allsamples2 %>% 
  filter(!grepl('ENSG', gene_stable_id))
nrow(x17321_allsamples2) # 17321(gene)*1017(samples) = 17615457
str(unique(x17321_allsamples2$gene_stable_id)) # 17237 gene
d<-as.data.table(table(x17321_allsamples2$gene_stable_id))
nrow(d)
head(d, 100)
table(d$N)
# 1017  2034  3051  4068 12204 
#17169    62     4     1     1
subset(d, N>=2034) 
x17321_dupli<-subset(d, N>=2034) 

'''
                 V1     N
 1: ENSG00000015568  3051
2: ENSG00000059145  2034
3: ENSG00000072310  2034
4: ENSG00000109181  2034
5: ENSG00000125551  2034
6: ENSG00000131548  2034
7: ENSG00000148483  2034
8: ENSG00000154864  4068
9: ENSG00000157600  2034
10: ENSG00000164871  2034
11: ENSG00000168477  2034
12: ENSG00000171489  2034
13: ENSG00000172058  2034
14: ENSG00000172288  2034
15: ENSG00000183385  2034
16: ENSG00000183675  2034
17: ENSG00000183753  3051
18: ENSG00000184110  2034
19: ENSG00000184324  2034
20: ENSG00000184750  2034
21: ENSG00000184779  2034
22: ENSG00000185565  2034
23: ENSG00000185684  2034
24: ENSG00000185751  2034
25: ENSG00000185945  2034
26: ENSG00000187243  2034
27: ENSG00000187627  2034
28: ENSG00000187653  2034
29: ENSG00000187775  3051
30: ENSG00000187786  2034
31: ENSG00000189013  2034
32: ENSG00000197681  2034
33: ENSG00000198040  2034
34: ENSG00000198393  2034
35: ENSG00000198457  2034
36: ENSG00000203811  2034
37: ENSG00000204164  2034
38: ENSG00000204344  2034
39: ENSG00000204828  2034
40: ENSG00000205176  2034
41: ENSG00000205944  2034
42: ENSG00000206258  2034
43: ENSG00000206338  2034
44: ENSG00000206342  2034
45: ENSG00000213599  2034
46: ENSG00000215765  2034
47: ENSG00000215790  2034
48: ENSG00000225932  2034
49: ENSG00000226033  2034
50: ENSG00000226257  2034
51: ENSG00000229341  2034
52: ENSG00000229353  2034
53: ENSG00000231608  2034
54: ENSG00000231852  2034
55: ENSG00000232414  2034
56: ENSG00000232948  2034
57: ENSG00000233151  2034
58: ENSG00000233323  2034
59: ENSG00000234947  2034
60: ENSG00000235134  2034
61: ENSG00000235173  2034
62: ENSG00000236221  2034
63: ENSG00000236236  2034
64: ENSG00000236250  2034
65: ENSG00000237763  3051
66: ENSG00000240428  2034
67: ENSG00000241476  2034
68: ENSG00000250589 12204
'''

######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
### check gene 16502 sampels

x16502<-subset(x, N>=16502 & N<16814)
nrow(x16502) #samples = 1807
head(x16502)
colnames(x16502) = c("submitted_sample_id", "gene_n")
head(x16502)
x16502_allsamples<-merge(x=x16502, y=data, by = "submitted_sample_id", all.x = TRUE)
nrow(x16502_allsamples) # 44182117
x16502_allsamples2<-x16502_allsamples[grep("ENSG", x16502_allsamples$gene_stable_id),]
x16502_allsamples2 %>% 
  filter(!grepl('ENSG', gene_stable_id))
nrow(x16502_allsamples2) # 30376281
str(unique(x16502_allsamples2$gene_stable_id)) # 16751 gene
d<-as.data.table(table(x16502_allsamples2$gene_stable_id))
nrow(d)
head(d, 100)
table(d$N)
#1721  1741  1748  1759  1765  1766  1769  1770  1771  1772  1774  1775  1776  1777 
#1     1     1     1     1     2     2     1     1     2     1     1     3     1 
#1778  1779  1780  1782  1783  1784  1785  1786  1787  1788  1789  1790  1791  1792 
#2     4     7     2     5     6     1     2     7     9     8     7     8     7 
#1793  1794  1795  1796  1797  1798  1799  1800  1801  1802  1803  1804  1805  1806 
#10    10    12    14    18    12    16    15    13    46    56    75   128   326 
#1807  3568  3610  3612  3614  5421  7228 
#15851     1     1     2    45     5     1
subset(d, N>=3614) 
'''
                 V1    N
 1: ENSG00000015568 5421
2: ENSG00000059145 3614
3: ENSG00000072310 3614
4: ENSG00000109181 3614
5: ENSG00000154864 7228
6: ENSG00000164871 3614
7: ENSG00000168477 3614
8: ENSG00000171489 3614
9: ENSG00000172058 3614
10: ENSG00000182330 5421
11: ENSG00000183558 3614
12: ENSG00000184779 3614
13: ENSG00000185565 3614
14: ENSG00000185684 3614
15: ENSG00000187243 3614
16: ENSG00000187653 3614
17: ENSG00000187775 5421
18: ENSG00000189013 3614
19: ENSG00000197932 5421
20: ENSG00000198040 3614
21: ENSG00000198393 3614
22: ENSG00000198457 3614
23: ENSG00000203811 3614
24: ENSG00000203989 3614
25: ENSG00000204164 3614
26: ENSG00000204344 3614
27: ENSG00000204828 3614
28: ENSG00000205944 3614
29: ENSG00000206258 3614
30: ENSG00000206338 3614
31: ENSG00000206342 3614
32: ENSG00000215790 3614
33: ENSG00000225932 3614
34: ENSG00000226033 3614
35: ENSG00000226257 3614
36: ENSG00000229341 3614
37: ENSG00000229353 3614
38: ENSG00000231608 3614
39: ENSG00000231852 3614
40: ENSG00000232414 3614
41: ENSG00000233151 3614
42: ENSG00000233323 3614
43: ENSG00000233630 5421
44: ENSG00000234947 3614
45: ENSG00000235134 3614
46: ENSG00000235173 3614
47: ENSG00000236221 3614
48: ENSG00000236236 3614
49: ENSG00000236250 3614
50: ENSG00000240428 3614
51: ENSG00000241476 3614
'''
x16502_dupli<-subset(d, N>=3614) 

####
x_dupli<-rbind(x11968_dupli, x16502_dupli, x16813_dupli, x17321_dupli)
head(x_dupli)
x_dupli2<-as.data.table(table(x_dupli$V1))
head(x_dupli2)
table(x_dupli2$N)
# 1  2  3  4 
#19  9 20 27
subset(x_dupli2, N==4) 
'''
                 V1 N
 1: ENSG00000015568 4
2: ENSG00000059145 4
3: ENSG00000072310 4
4: ENSG00000154864 4
5: ENSG00000168477 4
6: ENSG00000184779 4
7: ENSG00000185565 4
8: ENSG00000187775 4
9: ENSG00000189013 4
10: ENSG00000198040 4
11: ENSG00000198393 4
12: ENSG00000204344 4
13: ENSG00000206258 4
14: ENSG00000206342 4
15: ENSG00000215790 4
16: ENSG00000226033 4
17: ENSG00000226257 4
18: ENSG00000229341 4
19: ENSG00000229353 4
20: ENSG00000231608 4
21: ENSG00000233323 4
22: ENSG00000234947 4
23: ENSG00000235173 4
24: ENSG00000236221 4
25: ENSG00000236236 4
26: ENSG00000236250 4
27: ENSG00000241476 4
'''


########################################
# match gene id("ENSG") with gene name.#
########################################
#source("http://bioconductor.org/biocLite.R")
#biocLite("mygene")
install.packages('BiocManager')
BiocManager::install('mygene')



y<-as.data.table(table(sample_ENSG$gene_stable_id))
head(y)
summary(y$N)
table(y$N)
ENSG_list<-subset(y, N>=3948)

ENSG_11361<-subset(y, N==3948)
head(ENSG_11361)
nrow(ENSG_11361)

ENSG_27_dupli<-subset(x_dupli2, N==4) 
head(ENSG_27_dupli)

ENSG_11388<-rbind(ENSG_11361, ENSG_27_dupli)
nrow(ENSG_11388)

#1
library(mygene)
sample_ENSG_mach<-getGenes(ENSG_11388$V1, fields='symbol')
head(sample_ENSG_mach)
'''
             query  X_id  X_score   symbol notfound
1: ENSG00000000003  7105 18.66361   TSPAN6       NA
2: ENSG00000000005 64102 18.64555     TNMD       NA
3: ENSG00000000419  8813 18.71176     DPM1       NA
4: ENSG00000000457 57147 18.18913    SCYL3       NA
5: ENSG00000000460 55732 18.18781 C1orf112       NA
6: ENSG00000000938  2268 18.63222      FGR       NA
'''
sample_ENSG_mach<-as.data.table(sample_ENSG_mach)
sapply(sample_ENSG_mach, function(x) sum(is.na(x))) # count each columns 
#query     X_id  X_score   symbol notfound 
#   0      132      132      132    11256
nrow(sample_ENSG_mach)
sample_ENSG2<-sample_ENSG_mach[complete.cases(sample_ENSG_mach[ ,c("symbol")]), ]
sapply(sample_ENSG2, function(x) sum(is.na(x))) # count each columns 
head(sample_ENSG2)
nrow(sample_ENSG2)
str(unique(sample_ENSG2$symbol)) #10580
#######################################
#######################################
#######################################
sample_ENSG3<-as.data.table(table(sample_ENSG2$symbol))
nrow(sample_ENSG3)
head(sample_ENSG3)
table(sample_ENSG3$N)
#    1     2     3     4     5     6     7     8 
#10462     2     2     2    10    24    42    36

nrow(sample_ENSG3) #10580

head(sample_ENSG_na_genelist)
nrow(sample_ENSG_na_genelist) #18295

tumorsamples_ENSGna_ENSG_samegenelist<-(rbind(sample_ENSG_na_genelist, sample_ENSG3))
head(tumorsamples_ENSGna_ENSG_samegenelist)
table((as.data.table(table(tumorsamples_ENSGna_ENSG_samegenelist$V1)))$N)
####################################################################################33
###########################      normal samples    #################################33
####################################################################################33
library(RMySQL)
drv <- dbDriver("MySQL")
con <- dbConnect(MySQL(), 
                 dbname = "tcga", 
                 user = "solbi", 
                 host = "103.22.220.149", 
                 port = 3306,
                 password = "yonsei2018!")


r_sql <- c("SELECT `icgc_donor_id`,`project_code`,`submitted_sample_id`,`gene_stable_id`,`normalized_expression_level` FROM `gene_expression` WHERE `submitted_sample_id`like'%-11A%'")
data <- dbGetQuery(con, r_sql) #running 20min
data_SolidA<-data
r_sql <- c("SELECT `icgc_donor_id`,`project_code`,`submitted_sample_id`,`gene_stable_id`,`normalized_expression_level` FROM `gene_expression` WHERE `submitted_sample_id`like'%-11B%'")
data <- dbGetQuery(con, r_sql) #running 20min
data_SolidB<-data
r_sql <- c("SELECT `icgc_donor_id`,`project_code`,`submitted_sample_id`,`gene_stable_id`,`normalized_expression_level` FROM `gene_expression` WHERE `submitted_sample_id`like'%-11C%'")
data <- dbGetQuery(con, r_sql) #running 20min
data_SolidC<-data

data_Solid<-rbind(data_SolidA, data_SolidB, data_SolidC)
head(data_Solid)

head(table(data_Solid$submitted_sample_id), 10)
a<-as.data.table(table(data_Solid$submitted_sample_id))
summary(a$N)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#16796   18457   18457   19500   18457   35270
table(a$N)
#16796 16800 16806 16809 16810 16811 16812 16813 18457 35248 35255 35258 35259 35261 35262 35264 35265 35266 35267 35268 35269 35270 
#    1     1     1     1     1     1     6    11   515     1     1     1     1     1     2     2     1     3     4     1     7    13 
nrow(a) #576

##### exclude "ENSG"
library(dplyr)
normalsample_ENSG_na<-data_Solid %>% 
  filter(!grepl('ENSG', gene_stable_id))
head(normalsample_ENSG_na, 20)
normalsample_ENSG_na[grep("ENSG", normalsample_ENSG_na$gene_stable_id),]

str(unique(normalsample_ENSG_na$submitted_sample_id)) #553/576
str(unique(normalsample_ENSG_na$gene_stable_id)) #18456
nrow(normalsample_ENSG_na) #10206721
x<-as.data.table(table(normalsample_ENSG_na$submitted_sample_id))
summary(x$N)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#18457   18457   18457   18457   18457   18457

table(x$N)
#18457 
#  553

y<-as.data.table(table(normalsample_ENSG_na$gene_stable_id))
head(y)
summary(y$N)
nrow(y)
table(y$N) #
#  553  1106 
#18455     1 

subset(y, N==1106)
#        V1    N
#1: SLC35E2 1106

nrow(y)
normalsample_ENSG_na_genelist<-y #18456

#####  include "ENSG"
normalsample_ENSG<-data_Solid[grep("ENSG", data_Solid$gene_stable_id),]
normalsample_ENSG %>% 
  filter(!grepl('ENSG', gene_stable_id))

str(unique(normalsample_ENSG$submitted_sample_id)) #61/576
str(unique(normalsample_ENSG$gene_stable_id)) #16751

head(normalsample_ENSG)

x<-as.data.table(table(normalsample_ENSG$submitted_sample_id))
head(x)
summary(x$N)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#16791   16809   16812   16810   16813   16813
table(x$N)
#16791 16796 16798 16800 16801 16802 16804 16805 16806 16807 16808 16809 16810 16811 16812 16813 
#    1     1     1     1     1     1     1     2     1     2     1     4     5     2    13    24 

y<-as.data.table(table(normalsample_ENSG$gene_stable_id))
head(y)
summary(y$N)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#57.00   61.00   61.00   61.21   61.00  244.00 
table(y$N)
#57    58    59    60    61   122   183   244 
# 2     6    18   125 16545    49     5     1


x16502<-x
nrow(x16502) #samples = 61
head(x16502)
colnames(x16502) = c("submitted_sample_id", "gene_n")
head(x16502)
x16502_allsamples<-merge(x=x16502, y=data_Solid, by = "submitted_sample_id", all.x = TRUE)
nrow(x16502_allsamples) # 1726772
x16502_allsamples2<-x16502_allsamples[grep("ENSG", x16502_allsamples$gene_stable_id),]
x16502_allsamples2 %>% 
  filter(!grepl('ENSG', gene_stable_id))
nrow(x16502_allsamples2) # 1025406
str(unique(x16502_allsamples2$gene_stable_id)) # 16751 gene
d<-as.data.table(table(x16502_allsamples2$gene_stable_id))
nrow(d)
head(d)
table(d$N)
#57    58    59    60    61   122   183   244 
# 2     6    18   125 16545    49     5     1 
subset(d, N>=122)
'''
                 V1   N
 1: ENSG00000015568 183
2: ENSG00000059145 122
3: ENSG00000072310 122
4: ENSG00000109181 122
5: ENSG00000154864 244
6: ENSG00000164871 122
7: ENSG00000168477 122
8: ENSG00000171489 122
9: ENSG00000172058 122
10: ENSG00000172352 122
11: ENSG00000182330 183
12: ENSG00000183385 122
13: ENSG00000183558 122
14: ENSG00000184779 122
15: ENSG00000185565 122
16: ENSG00000185684 122
17: ENSG00000187243 122
18: ENSG00000187627 122
19: ENSG00000187653 122
20: ENSG00000187775 183
21: ENSG00000189013 122
22: ENSG00000197932 183
23: ENSG00000198040 122
24: ENSG00000198393 122
25: ENSG00000198457 122
26: ENSG00000203811 122
27: ENSG00000203989 122
28: ENSG00000204164 122
29: ENSG00000204344 122
30: ENSG00000204828 122
31: ENSG00000205176 122
32: ENSG00000205944 122
33: ENSG00000206258 122
34: ENSG00000206338 122
35: ENSG00000206342 122
36: ENSG00000215790 122
37: ENSG00000225932 122
38: ENSG00000226033 122
39: ENSG00000226257 122
40: ENSG00000229341 122
41: ENSG00000229353 122
42: ENSG00000231608 122
43: ENSG00000231852 122
44: ENSG00000232414 122
45: ENSG00000233151 122
46: ENSG00000233323 122
47: ENSG00000233630 183
48: ENSG00000234947 122
49: ENSG00000235134 122
50: ENSG00000235173 122
51: ENSG00000236221 122
52: ENSG00000236236 122
53: ENSG00000236250 122
54: ENSG00000240428 122
55: ENSG00000241476 122
'''
###
x16502_dupli<-subset(d, N>=3614) 

####
x_dupli<-subset(d, N>=122)
head(x_dupli)
nrow(x_dupli) #55
x_dupli2<-as.data.table(table(x_dupli$V1))

########################################
# match gene id("ENSG") with gene name.#
########################################
#source("http://bioconductor.org/biocLite.R")
#biocLite("mygene")
y<-as.data.table(table(normalsample_ENSG$gene_stable_id))
head(y)
summary(y$N)
table(y$N)
#57    58    59    60    61   122   183   244 
# 2     6    18   125 16545    49     5     1 

ENSG_16545<-subset(y, N==61)
head(ENSG_16545)
nrow(ENSG_16545)

ENSG_55_dupli<-x_dupli2
head(ENSG_55_dupli)
nrow(ENSG_55_dupli)
ENSG_11388<-rbind(ENSG_16545, ENSG_55_dupli)

#1
library(mygene)
normalsample_ENSG<-getGenes(ENSG_11388$V1, fields='symbol')
head(normalsample_ENSG)
'''
DataFrame with 6 rows and 5 columns
            query         _id   X_score      symbol  notfound
      <character> <character> <numeric> <character> <logical>
1 ENSG00000000003        7105  18.64988      TSPAN6        NA
2 ENSG00000000005       64102  17.78340        TNMD        NA
3 ENSG00000000419        8813  18.73387        DPM1        NA
4 ENSG00000000457       57147  18.18171       SCYL3        NA
5 ENSG00000000460       55732  18.18781    C1orf112        NA
6 ENSG00000000938        2268  17.48855         FGR        NA

'''
normalsample_ENSG<-as.data.table(normalsample_ENSG)
sapply(normalsample_ENSG, function(x) sum(is.na(x))) # count each columns 
#query     X_id  X_score   symbol notfound 
#0      256      256      256    16344 


#######################################
#######################################
#######################################
normalsample_ENSG3<-as.data.table(table(normalsample_ENSG$symbol))
nrow(normalsample_ENSG3)
head(normalsample_ENSG3)
table(normalsample_ENSG3$N)
#    1     2     3     4     5     6     7     8 
#15334     4     2     3    11    28    55    47

nrow(normalsample_ENSG3) #15484

head(normalsample_ENSG_na_genelist)
nrow(normalsample_ENSG_na_genelist) #18456

normalsamples_ENSGna_ENSG_samegenelist<-(rbind(normalsample_ENSG_na_genelist, sample_ENSG3))
head(normalsamples_ENSGna_ENSG_samegenelist)
table((as.data.table(table(normalsamples_ENSGna_ENSG_samegenelist$V1)))$N)


################################################################################################
##########################   make dataset use machine learninig ################################
################################################################################################

########  1. primary dataset

##############################
###check gene 18457 sampels
x18457<-subset(x, N==18457)
head(x18457)
colnames(x18457) = c("submitted_sample_id", "gene_n")
head(x18457)
x18457_allsamples<-merge(x=x18457, y=data, by = "submitted_sample_id", all.x = TRUE)
nrow(x18457_allsamples) # 113939354
x18457_allsamples2<-x18457_allsamples %>% 
  filter(!grepl('ENSG', gene_stable_id))
x18457_allsamples2[grep("ENSG", x18457_allsamples2$gene_stable_id),]

nrow(x18457_allsamples2) # 18457(gene)*4838(samples) = 101365844
str(unique(x18457_allsamples2$gene_stable_id)) #18456 ("SLC35E2" expression data is 2)
d<-as.data.table(table(x18457_allsamples2$gene_stable_id))
head(d)
table(d$N)
subset(d, N==10984) 

#        V1     N
#1: SLC35E2 10984
''

head(x18457_allsamples2)  
nrow(x18457_allsamples2) # 18457(gene)*4838(samples)=101365844
x18457_allsamples_SLC35E2<-subset(x18457_allsamples2, gene_stable_id=="SLC35E2")

nrow(x18457_allsamples_SLC35E2) # 4838(samples)*2=10984
head(x18457_allsamples_SLC35E2) #

write.csv(x18457_allsamples2, file = "/home/hwanghou/20200311/tcga_primary_4838sample_18456gene_rowdata.csv", row.names = FALSE)

######################################################################################
################make table x=submitted_sample_id, y=gene_stable_id####################
######################################################################################

#install.packages("reshape2")
tcga_primary_dataset_4838<-read.csv("/home/hwanghou/20200311/tcga_primary_4838sample_18456gene_rowdata.csv")

library(reshape2)
head(tcga_primary_dataset_4838)
tcga_primary_dataset_4838_geneexpres<-tcga_primary_dataset_4838[c("submitted_sample_id", "gene_stable_id", "normalized_expression_level")]
head(tcga_primary_dataset_4838_geneexpres)
str(unique(tcga_primary_dataset_4838_geneexpres$submitted_sample_id)) #4838


tcga_primary_dataset_4838_geneexpres_dataset<-dcast(tcga_primary_dataset_4838_geneexpres, submitted_sample_id~gene_stable_id, value.var = "normalized_expression_level", sum)
ncol(tcga_primary_dataset_4838_geneexpres_dataset) #18457
nrow(tcga_primary_dataset_4838_geneexpres_dataset) #4838

#samples feature data 
tcgabatabase_feature<-read.csv("/home/hwanghou/20200311/tcga_patient (5).csv")
head(tcgabatabase_feature)
nrow(tcgabatabase_feature)
tcgabatabase_feature2<-unique(tcgabatabase_feature)
nrow(tcgabatabase_feature2)
table(tcgabatabase_feature2$pathologic_stage)
table(tcgabatabase_feature2$clinical_stage)
table(tcgabatabase_feature2$race)
tcgabatabase_feature2 <- tcgabatabase_feature2[,c(2,5)]
names(tcgabatabase_feature2)[1] <- c('submitted_sample_id')

#samples clinical data
tcgabatabase_clinical <- read.csv("/home/hwanghou/20200311/clinical (2).csv")
head(tcgabatabase_clinical)
tcgabatabase_clinical <- tcgabatabase_clinical[,c(1,2,4,8)]
names(tcgabatabase_clinical)[1]<- c('submitted_sample_id')
names(tcgabatabase_clinical)[2]<- c('age')

tcga_primary_dataset_4838_geneexpres_dataset$submitted_sample_id<-substr(tcga_primary_dataset_4838_geneexpres_dataset$submitted_sample_id,1,12)
names(tcgabatabase_feature2)[2] <- c('submitted_sample_id')
head(tcga_primary_dataset_4838_geneexpres_dataset)


stage <- read.table('/home/hwanghou/20200312/Final_clinical.tsv',header=TRUE,sep='\t')
table(stage$Stage)
table(stage$Gender)
head(stage)
stage <- stage[,c(2,4,5,7,9)]
nrow(stage)
head(stage)
stage<-unique(stage)
names(stage)[1] <- c('submitted_sample_id')

tcgadatabase_patient_feature <- read.csv('/home/hwanghou/20200313/tcgadatabase_patient_feature.csv')
tcgadatabase_patient_feature <- read.csv('/home/hwanghou/20200316/tcgadatabase_patient_feature_3.csv')
names(tcgadatabase_patient_feature)[1] <- c('submitted_sample_id')
head(tcgadatabase_patient_feature)



tcga_primary_dataset_4838_geneexpres_dataset2 <- merge(x = tcgadatabase_patient_feature, y = tcga_primary_dataset_4838_geneexpres_dataset, by = 'submitted_sample_id', all.y = TRUE)
tcga_primary_dataset_4838_geneexpres_dataset2 <- merge(x=tcgabatabase_feature2,y=tcga_primary_dataset_4838_geneexpres_dataset2,by = 'submitted_sample_id',all.y=TRUE)
tcga_primary_dataset_4838_geneexpres_dataset2_unique <- unique(tcga_primary_dataset_4838_geneexpres_dataset2)
tcga_primary_dataset_4838_geneexpres_dataset2_unique <- tcga_primary_dataset_4838_geneexpres_dataset2_unique[,-c(5)]


tcga_primary_dataset_4838_geneexpres_dataset2 <- merge(x = tcgabatabase_feature2, y = tcga_primary_dataset_4838_geneexpres_dataset, by = 'submitted_sample_id', all.y = TRUE)
tcga_primary_dataset_4838_geneexpres_dataset2 <- merge(x=tcgabatabase_clinical,y=tcga_primary_dataset_4838_geneexpres_dataset2,by = 'submitted_sample_id',all.y=TRUE)
tcga_primary_dataset_4838_geneexpres_dataset2_unique <- unique(tcga_primary_dataset_4838_geneexpres_dataset2)
tcga_primary_dataset_4838_geneexpres_dataset2_unique <- tcga_primary_dataset_4838_geneexpres_dataset2_unique[,-c(5)]

nrow(tcga_primary_dataset_4838_geneexpres_dataset)
ncol(tcga_primary_dataset_4838_geneexpres_dataset)
nrow(tcga_primary_dataset_4838_geneexpres_dataset2)
ncol(tcga_primary_dataset_4838_geneexpres_dataset2)
head(colnames(tcga_primary_dataset_4838_geneexpres_dataset2_unique), 30)
table(tcga_primary_dataset_4838_geneexpres_dataset2_unique$Race)


write.csv(tcga_primary_dataset_4838_geneexpres_dataset2_unique, file = "/home/hwanghou/20200311/tcga_primary_dataset_4838_geneexpres_dataset.csv", row.names = FALSE)
tcga_primary_dataset_4838_geneexpres_dataset2<-read.csv("/home/hwanghou/20200311/tcga_primary_dataset_4838_geneexpres_dataset.csv")
sapply(tcga_primary_dataset_4838_geneexpres_dataset2_unique, function(x) sum(is.na(x))) # count each columns 


'''
submitted_sample_id                 icgc_donor_id            submitted_donor_id                  project_code                     donor_sex        donor_age_at_diagnosis 
0                             0                             0                             0                             0                             0 
bcr_patient_barcode                          race           project_code_rename              donor_sex_rename donor_age_at_diagnosis_rename                   race_rename 
0                             1                             0                             0                             0                             1 
Status_rename                  Stage_rename                          A1BG                          A1CF                           A2M                         A2ML1 
1516                          1553
'''

nrow(tcga_primary_dataset_4838_geneexpres_dataset)
nrow(tcga_primary_dataset_4838_geneexpres_dataset2)
table(tcga_primary_dataset_4838_geneexpres_dataset2_unique$donor_age_at_diagnosis_rename)
table(tcga_primary_dataset_4838_geneexpres_dataset2_unique$donor_sex_rename)
tcga_primary_dataset_4838_geneexpres_dataset2_unique$donor_sex_rename <- gsub('male',1,tcga_primary_dataset_4838_geneexpres_dataset2_unique$donor_sex_rename)
#    0    1 
# 3242 2250 
table(tcga_primary_dataset_4838_geneexpres_dataset2_unique$race_rename)
#[Not Available] [Not Evaluated]       [Unknown]               0               1               2               3               4 
#            577              36              39              13             242             416               9            4159 

table(tcga_primary_dataset_4838_geneexpres_dataset2_unique$Stage_rename)
table(tcga_primary_dataset_4838_geneexpres_dataset2$pathologic_stage)
#   0         1         2         3 Stage Tis   Stage X        T2        T3        T4 
#1381      1072      1021       308         1         8        65        78         5
table(tcga_primary_dataset_4838_geneexpres_dataset2_unique$donor_age_at_diagnosis)
# 0    1    2    3    4    5    6    7    8 
#51  114  336  716 1272 1553 1090  337   23 
table(tcga_primary_dataset_4838_geneexpres_dataset2$donor_age_at_diagnosis_rename)

table(tcga_primary_dataset_4838_geneexpres_dataset2_unique$project_code_rename)
#   0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  18  19 
# 182 938  61 376 142 351 493 126 255 121 442 396 252  50 173 140  57 474 463 
table(tcga_primary_dataset_4838_geneexpres_dataset2$project_code_rename)

str(unique(tcga_primary_dataset_4838_geneexpres_dataset2$race_rename))
sum(is.na(tcga_primary_dataset_4838_geneexpres_dataset2_unique$Stage_rename))
sum(is.na(tcga_primary_dataset_4838_geneexpres_dataset2_unique$race_rename))
sum(is.na(tcga_primary_dataset_4838_geneexpres_dataset2_unique$donor_age_at_diagnosis))

str(tcga_primary_dataset_4838_geneexpres_dataset2)
tcga_primary_dataset_4838_geneexpres_dataset2 <- tcga_primary_dataset_4838_geneexpres_dataset2_unique


tcga_primary_dataset_4838_geneexpres_dataset2$project_code_rename<-gsub("18", "17", tcga_primary_dataset_4838_geneexpres_dataset2$project_code_rename)
tcga_primary_dataset_4838_geneexpres_dataset2$project_code_rename<-gsub("19", "18", tcga_primary_dataset_4838_geneexpres_dataset2$project_code_rename)
tcga_primary_dataset_4838_geneexpres_dataset2 <- droplevels(tcga_primary_dataset_4838_geneexpres_dataset2)
table(tcga_primary_dataset_4838_geneexpres_dataset2$project_code_rename)
write.csv(tcga_primary_dataset_4838_geneexpres_dataset4, file = "/home/hwanghou/20200313//tcga_primary_dataset_4838_geneexpres_dataset.csv", row.names = FALSE)


###### first table 
table(tcga_primary_dataset_4838_geneexpres_dataset2$donor_age_at_diagnosis_rename) 
table(tcga_primary_dataset_4838_geneexpres_dataset2$donor_age_at_diagnosis) 
'''
0  14  15  16  17  18  19  20  21  22  23  24  25  26  27  28  29  30  31  32  33  34  35  36  37  38  39  40  41  42  43  44  45  46  47  48 
37   1   1   1   2   4   5   8   7   4  12  12   7  17   9  19  19  27  25  20  29  41  28  34  39  50  43  63  58  54  63  44  84  78  95  85 
49  50  51  52  53  54  55  56  57  58  59  60  61  62  63  64  65  66  67  68  69  70  71  72  73  74  75  76  77  78  79  80  81  82  83  84 
92  91 126 111 129 124 113 132 148 149 149 180 161 177 167 141 151 145 146 144 141 133 127 124 131 121 111  98  93  75  77  70  50  46  37  52 
85  86  87  88  89  90 
27  14  19  14   8  23
'''
tcga_primary_dataset_4838_geneexpres_dataset3<-tcga_primary_dataset_4838_geneexpres_dataset2[tcga_primary_dataset_4838_geneexpres_dataset2$donor_age_at_diagnosis!=0, ]
tcga_primary_dataset_4838_geneexpres_dataset2[tcga_primary_dataset_4838_geneexpres_dataset2$donor_age_at_diagnosis==0,c(0:6)]

table(tcga_primary_dataset_4838_geneexpres_dataset3$donor_age_at_diagnosis) 
table(tcga_primary_dataset_4838_geneexpres_dataset3$donor_age_at_diagnosis_rename) 
# 0    1    2    3    4    5    6    7    8 
#14  114  336  716 1272 1553 1090  337   23   

nrow(tcga_primary_dataset_4838_geneexpres_dataset3)
tcga_primary_dataset_4838_geneexpres_dataset4<-tcga_primary_dataset_4838_geneexpres_dataset3[complete.cases(tcga_primary_dataset_4838_geneexpres_dataset3[ ,c("race_rename")]), ]# remove na = 16 samples
nrow(tcga_primary_dataset_4838_geneexpres_dataset4)
table(tcga_primary_dataset_4838_geneexpres_dataset4$race_rename)
#[Not Available] [Not Evaluated]       [Unknown]               0               1               2               3               4 
#           547              30              39              13             242             416               9            4158  

tcga_primary_dataset_4838_geneexpres_dataset4<-tcga_primary_dataset_4838_geneexpres_dataset4[tcga_primary_dataset_4838_geneexpres_dataset4$race_rename!="[Not Evaluated]" & tcga_primary_dataset_4838_geneexpres_dataset4$race_rename!="[Not Available]" & tcga_primary_dataset_4838_geneexpres_dataset4$race_rename!="[Unknown]", ]
table(tcga_primary_dataset_4838_geneexpres_dataset4$race_rename)
tcga_primary_dataset_4838_geneexpres_dataset4 <- droplevels(tcga_primary_dataset_4838_geneexpres_dataset4)
table(tcga_primary_dataset_4838_geneexpres_dataset4$race_rename)
# 0    1    2    3    4 
#13  242  416    9 4158

table(tcga_primary_dataset_4838_geneexpres_dataset4$project_code_rename)
str(unique(tcga_primary_dataset_4838_geneexpres_dataset4$project_code_rename))
table(tcga_primary_dataset_4838_geneexpres_dataset4$project_code)
table(tcga_primary_dataset_4838_geneexpres_dataset4$donor_sex_rename)

nrow(tcga_primary_dataset_4838_geneexpres_dataset4)

write.csv(tcga_primary_dataset_4838_geneexpres_dataset4, file = "/home/hwanghou/20200316/tcga_primary_dataset_4838_geneexpres_dataset_4featuredata.csv", row.names = FALSE)

table(tcga_primary_dataset_4838_geneexpres_dataset4$Stage_rename)
#   0         1         2         3 Stage Tis   Stage X        T2        T3        T4 
#1381      1072      1021       308         1         8        65        78         5 

table(tcga_primary_dataset_4838_geneexpres_dataset4$project_code_rename)

sapply(tcga_primary_dataset_4838_geneexpres_dataset4, function(x) sum(is.na(x))) # count each columns
'''
submitted_sample_id                 icgc_donor_id            submitted_donor_id                  project_code                     donor_sex 
0                             0                             0                             0                             0 
donor_age_at_diagnosis           bcr_patient_barcode                          race           project_code_rename              donor_sex_rename 
0                             0                             0                             0                             0 
donor_age_at_diagnosis_rename                   race_rename                 Status_rename                  Stage_rename                          A1BG 
0                             0                           862                           899                             0 
'''

tcga_primary_dataset_4838_geneexpres_dataset5<-tcga_primary_dataset_4838_geneexpres_dataset4[complete.cases(tcga_primary_dataset_4838_geneexpres_dataset4[ ,c("Stage_rename")]), ]# remove na = 899
nrow(tcga_primary_dataset_4838_geneexpres_dataset4)
nrow(tcga_primary_dataset_4838_geneexpres_dataset5)
table(tcga_primary_dataset_4838_geneexpres_dataset5$Stage_rename)
#   0         1         2         3 Stage Tis   Stage X        T2        T3        T4 
#1381      1072      1021       308         1         8        65        78         5 
tcga_primary_dataset_4838_geneexpres_dataset5<-tcga_primary_dataset_4838_geneexpres_dataset5[tcga_primary_dataset_4838_geneexpres_dataset5$Stage_rename!="Stage Tis" & tcga_primary_dataset_4838_geneexpres_dataset5$Stage_rename!="Stage X", ]
table(tcga_primary_dataset_4838_geneexpres_dataset5$Stage_rename)
#   0         1         2         3 Stage Tis   Stage X        T2        T3        T4 
#1381      1072      1021       308         0         0         0         0         0
tcga_primary_dataset_4838_geneexpres_dataset5 <- droplevels(tcga_primary_dataset_4838_geneexpres_dataset5)
table(tcga_primary_dataset_4838_geneexpres_dataset5$Stage_rename)

#   0    1    2    3 
#1381 1072 1021  308 
nrow(tcga_primary_dataset_4838_geneexpres_dataset5)
head(colnames(tcga_primary_dataset_4838_geneexpres_dataset4), 30)

tcga_primary_dataset_4838_geneexpres_dataset5$donor_sex_rename <- gsub('male',1,tcga_primary_dataset_4838_geneexpres_dataset5$donor_sex_rename)
write.csv(tcga_primary_dataset_4838_geneexpres_dataset5, file = "/home/hwanghou/20200316/tcga_primary_dataset_3782_geneexpres_dataset_stage_featuredata.csv", row.names = FALSE)


##### normal dataset
library(RMySQL)
drv <- dbDriver("MySQL")
con <- dbConnect(MySQL(), 
                 dbname = "tcga", 
                 user = "solbi", 
                 host = "103.22.220.149", 
                 port = 3306,
                 password = "yonsei2018!")


r_sql <- c("SELECT `icgc_donor_id`,`project_code`,`submitted_sample_id`,`gene_stable_id`,`normalized_expression_level` FROM `gene_expression` WHERE `submitted_sample_id`like'%-11A%'")
data <- dbGetQuery(con, r_sql) #running 20min
data_SolidA<-data
r_sql <- c("SELECT `icgc_donor_id`,`project_code`,`submitted_sample_id`,`gene_stable_id`,`normalized_expression_level` FROM `gene_expression` WHERE `submitted_sample_id`like'%-11B%'")
data <- dbGetQuery(con, r_sql) #running 20min
data_SolidB<-data
r_sql <- c("SELECT `icgc_donor_id`,`project_code`,`submitted_sample_id`,`gene_stable_id`,`normalized_expression_level` FROM `gene_expression` WHERE `submitted_sample_id`like'%-11C%'")
data <- dbGetQuery(con, r_sql) #running 20min
data_SolidC<-data

data_Solid<-rbind(data_SolidA, data_SolidB, data_SolidC)
head(data_Solid)

head(table(data_Solid$submitted_sample_id), 10)
a<-as.data.table(table(data_Solid$submitted_sample_id))
summary(a$N)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#16796   18457   18457   19500   18457   35270
table(a$N)
#16796 16800 16806 16809 16810 16811 16812 16813 18457 35248 35255 35258 35259 35261 35262 35264 35265 35266 35267 35268 35269 35270 
#    1     1     1     1     1     1     6    11   515     1     1     1     1     1     2     2     1     3     4     1     7    13 
nrow(a) #576

##### exclude "ENSG"
library(dplyr)
normalsample_ENSG_na<-data_Solid %>% 
  filter(!grepl('ENSG', gene_stable_id))
head(normalsample_ENSG_na, 20)
normalsample_ENSG_na[grep("ENSG", normalsample_ENSG_na$gene_stable_id),]

str(unique(normalsample_ENSG_na$submitted_sample_id)) #553/576
str(unique(normalsample_ENSG_na$gene_stable_id)) #18456
nrow(normalsample_ENSG_na) #10206721
x<-as.data.table(table(normalsample_ENSG_na$submitted_sample_id))
summary(x$N)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#18457   18457   18457   18457   18457   18457

table(x$N)
#18457 
#  553
write.csv(normalsample_ENSG_na, file = "/home/hwanghou/20200311/tcga_normal_553sample_18456gene_rowdata.csv", row.names = FALSE)

######################################################################################
################make table x=submitted_sample_id, y=gene_stable_id####################
######################################################################################

#install.packages("reshape2")
tcga_normal_dataset_553 <- normalsample_ENSG_na
#tcga_normal_dataset_553<-read.csv("/home/solbi/20190118/tcga_normal_553sample_18456gene_rowdata.csv")
library(reshape2)
head(tcga_normal_dataset_553)
tcga_normal_dataset_553_geneexpres<-tcga_normal_dataset_553[c("submitted_sample_id", "gene_stable_id", "normalized_expression_level")]
head(tcga_normal_dataset_553_geneexpres)
str(unique(tcga_normal_dataset_553_geneexpres$submitted_sample_id)) #553


tcga_normal_dataset_553_geneexpres_dataset<-dcast(tcga_normal_dataset_553_geneexpres, submitted_sample_id~gene_stable_id, value.var = "normalized_expression_level", sum)
ncol(tcga_normal_dataset_553_geneexpres_dataset) #18457
nrow(tcga_normal_dataset_553_geneexpres_dataset) #553
write.csv(tcga_normal_dataset_553_geneexpres_dataset,file = "/home/hwanghou/20200312/tcga_normal_dataset_553_geneexpres_dataset.csv",row.names = FALSE)

#samples feature data 
tcgabatabase_feature<-read.csv("/home/solbi/20190118/tcgadatabase_patient_feature.csv")
head(tcgabatabase_feature)
nrow(tcgabatabase_feature)
tcgabatabase_feature2<-unique(tcgabatabase_feature)
nrow(unique(tcgabatabase_feature2$icgc_donor_id))
names(tcgabatabase_feature2[6]) <- c('submitted_sample_id')


#sample_feature
tcgadatabase_patient_feature <- read.csv('/home/hwanghou/20200316/tcgadatabase_patient_feature_3.csv')
names(tcgadatabase_patient_feature)[1] <- c('submitted_sample_id')
head(tcgadatabase_patient_feature)

tcga_normal_dataset_553_geneexpres_dataset$submitted_sample_id <- substr(tcga_normal_dataset_553_geneexpres_dataset2$submitted_sample_id,1,12)

tcga_normal_dataset_553_geneexpres_dataset2 <- merge(x = tcgadatabase_patient_feature, y = tcga_normal_dataset_553_geneexpres_dataset, by = 'submitted_sample_id', all.y = TRUE)
nrow(tcga_normal_dataset_553_geneexpres_dataset)
ncol(tcga_normal_dataset_553_geneexpres_dataset)
nrow(tcga_normal_dataset_553_geneexpres_dataset2)
ncol(tcga_normal_dataset_553_geneexpres_dataset2)
head(colnames(tcga_normal_dataset_553_geneexpres_dataset2), 30)

write.csv(tcga_normal_dataset_553_geneexpres_dataset2, file = "/home/hwnaghou/20200316/tcga_normal_dataset_553_geneexpres_dataset.csv", row.names = FALSE)
#tcga_normal_dataset_553_geneexpres_dataset2<-read.csv("/home/solbi/20190118/tcga_normal_dataset_553_geneexpres_dataset.csv")
sapply(tcga_normal_dataset_553_geneexpres_dataset2, function(x) sum(is.na(x))) # count each columns 

'''
          submitted_sample_id                 icgc_donor_id            submitted_donor_id                  project_code                     donor_sex 
                            0                             0                             0                             0                             0 
       donor_age_at_diagnosis           bcr_patient_barcode                          race           project_code_rename              donor_sex_rename 
                            0                             0                             1                             0                             0 
donor_age_at_diagnosis_rename                   race_rename                 Status_rename                  Stage_rename                          A1BG 
                            0                             1                           553                           553 
'''

nrow(tcga_normal_dataset_553_geneexpres_dataset)
nrow(tcga_normal_dataset_553_geneexpres_dataset2)
tcga_normal_dataset_553_geneexpres_dataset2 <- unique(tcga_normal_dataset_553_geneexpres_dataset2)
tcga_normal_dataset_553_geneexpres_dataset2$donor_sex_rename <- gsub('male',1,tcga_normal_dataset_553_geneexpres_dataset2$donor_sex_rename)
table(tcga_normal_dataset_553_geneexpres_dataset2$donor_sex_rename)

#    0    1 
#  293  260 
table(tcga_normal_dataset_553_geneexpres_dataset2$race_rename)
#[Not Available] [Not Evaluated]       [Unknown]               0               1               2               3               4 
#            38               1               1               1              11              41               0             459 

table(tcga_normal_dataset_553_geneexpres_dataset2$donor_age_at_diagnosis)
#  0   1   2   3   4   5   6   7   8 
#  3  12  31  74 108 162 116  42   5

table(tcga_normal_dataset_553_geneexpres_dataset2$donor_age_at_diagnosis_rename)

table(tcga_normal_dataset_553_geneexpres_dataset2$project_code_rename)
# 0   1   2   3   5   6   7   9  10  11  13  14  15  18  19 
#16 108   3  23  39  72  28  46  55  44   1  38   5  53  22 
table(tcga_normal_dataset_553_geneexpres_dataset2$project_code)
#BLCA-US BOCA-UK BRCA-UK BRCA-US CESC-US CLLE-ES CMDI-UK COAD-US EOPC-DE ESAD-UK  GBM-US HNSC-US KIRC-US KIRP-US LAML-US  LGG-US LICA-FR LIHC-US LINC-JP LIRI-JP LUAD-US 
#     16       0       0     108       3       0       0      23       0       0       0      39      72      28       0       0       0      46       0       0      55 
#LUSC-US MALY-DE  NBL-US ORCA-IN   OV-AU   OV-US PAAD-US PACA-AU PACA-CA PAEN-AU PBCA-DE PRAD-CA PRAD-US READ-US RECA-CN RECA-EU SKCM-US STAD-US THCA-SA THCA-US UCEC-US 
#     44       0       0       0       0       0       1       0       0       0       0       0      38       5       0       0       0       0       0      53      22 

str(table(tcga_normal_dataset_553_geneexpres_dataset2$project_code_rename)) #15

tcga_normal_dataset_553_geneexpres_dataset3<-tcga_normal_dataset_553_geneexpres_dataset2
table(tcga_normal_dataset_553_geneexpres_dataset2$project_code_rename)


nrow(tcga_normal_dataset_553_geneexpres_dataset2)
write.csv(tcga_normal_dataset_553_geneexpres_dataset2, file = "/home/hwanghou/20200316/tcga_normal_dataset_553_geneexpres_dataset.csv", row.names = FALSE)


###### first table 
table(tcga_normal_dataset_553_geneexpres_dataset2$project_code_rename)
tcga_normal_dataset_553_geneexpres_dataset3<-tcga_normal_dataset_553_geneexpres_dataset2
tcga_normal_dataset_553_geneexpres_dataset3$project_code_rename<-gsub("18", "17", tcga_normal_dataset_553_geneexpres_dataset3$project_code_rename)
tcga_normal_dataset_553_geneexpres_dataset3$project_code_rename<-gsub("19", "18", tcga_normal_dataset_553_geneexpres_dataset3$project_code_rename)
tcga_normal_dataset_553_geneexpres_dataset3 <- droplevels(tcga_normal_dataset_553_geneexpres_dataset3)
table(tcga_normal_dataset_553_geneexpres_dataset3$project_code_rename)

table(tcga_normal_dataset_553_geneexpres_dataset3$donor_age_at_diagnosis_rename) 
table(tcga_normal_dataset_553_geneexpres_dataset3$donor_age_at_diagnosis) 

tcga_normal_dataset_553_geneexpres_dataset4<-tcga_normal_dataset_553_geneexpres_dataset3
nrow(tcga_normal_dataset_553_geneexpres_dataset4)

tcga_normal_dataset_553_geneexpres_dataset4<-tcga_normal_dataset_553_geneexpres_dataset4[complete.cases(tcga_normal_dataset_553_geneexpres_dataset4[ ,c("race_rename")]), ]# remove na = 16 samples
nrow(tcga_normal_dataset_553_geneexpres_dataset4)
table(tcga_normal_dataset_553_geneexpres_dataset4$race_rename)
#[Not Available] [Not Evaluated]       [Unknown]               0               1               2               3 
#             38               1               1               1              11              41             459 

tcga_normal_dataset_553_geneexpres_dataset5<-tcga_normal_dataset_553_geneexpres_dataset4[tcga_normal_dataset_553_geneexpres_dataset4$race_rename!="[Not Evaluated]" & tcga_normal_dataset_553_geneexpres_dataset4$race_rename!="[Not Available]" & tcga_normal_dataset_553_geneexpres_dataset4$race_rename!="[Unknown]", ]
table(tcga_normal_dataset_553_geneexpres_dataset5$race_rename)
tcga_normal_dataset_553_geneexpres_dataset5 <- droplevels(tcga_normal_dataset_553_geneexpres_dataset5)
table(tcga_normal_dataset_553_geneexpres_dataset5$race_rename)
# 0    1    2    3    4 
#13  242  416    9 4158
table(tcga_normal_dataset_553_geneexpres_dataset5$project_code_rename)

nrow(tcga_normal_dataset_553_geneexpres_dataset5)
write.csv(tcga_normal_dataset_553_geneexpres_dataset5, file = "/home/hwanghou/20200316/tcga_normal_dataset_512_geneexpres_dataset_4featuredata.csv", row.names = FALSE)



###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################
####################################################          anova table            ##################################################################
###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################


tcga_primary_dataset_4838_geneexpres_dataset4<-read.csv("/home/hwanghou/20200312/tcga_primary_dataset_5825_geneexpres_dataset_4featuredata.csv")
tcga_primary_dataset_4838_geneexpres_dataset4
'''
[1] "submitted_sample_id"           "icgc_donor_id"                
[3] "submitted_donor_id"            "project_code"                 
[5] "donor_sex"                     "donor_age_at_diagnosis"       
[7] "bcr_patient_barcode"           "race"                         
[9] "project_code_rename"           "donor_sex_rename"             
[11] "donor_age_at_diagnosis_rename" "race_rename"                  
[13] "Status_rename"                 "Stage_rename"                 
[15] "A1BG"                          "A1CF"                         
[17] "A2M"                           "A2ML1"                        
[19] "A4GALT"                        "A4GNT"                        
[21] "AAAS"                          "AACS"                         
[23] "AADAC"                         "AADACL2"                      
[25] "AADACL3"   
'''
############################################################################################################################################################################
############################################################################################################################################################################
#########################################################          tcga all data t-test about gender           #############################################################
############################################################################################################################################################################
############################################################################################################################################################################
tcga<-read.csv("/home/solbi/20190118/tcga_primary_dataset_4838_geneexpres_dataset_4featuredata.csv")
tcga <- tcga_primary_dataset_4838_geneexpres_dataset4
#each genes t-test
library(data.table)
library(dplyr)
tcga_all_gender<-tcga[,-c(1, 2, 3, 4, 5,6, 7,8,9,11,12)]
head(colnames(tcga_all_gender))
#1
ttest_gender<-t(sapply(tcga_all_gender[-1], function(x) 
  +   unlist(anova(lm(x~donor_sex_rename, data=tcga_all_gender)))))
library(data.table)
ttest_gender<-as.data.table(ttest_gender, keep.rownames=TRUE) #library(data.table)
ttest_gender_bf<-cbind(ttest_gender, P.Bonferroni = p.adjust(ttest_gender$`Pr(>F)1`, method="bonferroni"))
write.csv(ttest_gender_bf, file = "/home/hwanghou/20200316/tcga_4838_ttest_gender.csv", row.names = FALSE)

                       
ttest_gender_bf<-read.csv("/home/solbi/20190118/tcga_4838_ttest_gender.csv")
head(ttest_gender_bf)
#table(ttest_gender_bf$P.Bonferroni)
#p.B = p.adjust(ttest_gender$p.value, method="bonferroni")
str(subset(ttest_gender_bf, P.Bonferroni<=0.05)) # 7757
str(subset(ttest_gender_bf, P.Bonferroni<=0.01)) # gene = 7220
str(subset(ttest_gender_bf, P.Bonferroni<=0.005)) # gene = 6998
str(subset(ttest_gender_bf, P.Bonferroni<=0.001)) # gene = 6508

###  p<0.05
gender<-tcga[,-c(1, 2, 3, 4, 5,6, 7,8,9,11,12)]
head(colnames(gender), 30)

table(gender$donor_sex_rename)

ttest_gender_bf<-read.csv("/home/solbi/20190118/tcga_4838_ttest_gender.csv")
head(ttest_gender_bf)
sapply(ttest_gender_bf, function(x) sum(is.na(x))) # count each columns 


anova_gender_p.B.05<-subset(ttest_gender_bf, P.Bonferroni<=0.05) # gene = 766
anova_gender_p.B.05 <- droplevels(anova_gender_p.B.05)
head(anova_gender_p.B.05)
str(anova_gender_p.B.05)

anova_gender_gene_list_05<-subset(anova_gender_p.B.05, select=c(rn))
colnames(anova_gender_gene_list_05) = c("gene_stable_id")
anova_gender_gene_list_05<-cbind(anova_gender_gene_list_05, donor_sex_rename=rep(1, nrow(anova_gender_gene_list_05)))
anova_gender_gene_list_05<-cbind(anova_gender_gene_list_05, num2=rep(2, nrow(anova_gender_gene_list_05)))
str(anova_gender_gene_list_05)
anova_gender_gene_list_05[, 2] <- as.numeric(as.character( anova_gender_gene_list_05[, 2] ))
anova_gender_gene_list_05[, 3] <- as.numeric(as.character( anova_gender_gene_list_05[, 3] ))
str(anova_gender_gene_list_05)

library(reshape2)
anova_gender_gene_list_table_05<-dcast(anova_gender_gene_list_05, donor_sex_rename ~ gene_stable_id, value.var = "num2", sum) #library(reshape2)
ncol(anova_gender_gene_list_table_05)
nrow(anova_gender_gene_list_table_05)
#write.csv(anova_gender_gene_list_table_05, file = "/home/solbi/20190118/4838anova_gender_gene_list_table_05.csv", row.names = FALSE)

#anova_gender_gene_list_table_05<-read.csv("/home/solbi/20190118/4838anova_gender_gene_list_table_05.csv")
str(anova_gender_gene_list_table_05)
samecols05 <- intersect(colnames(anova_gender_gene_list_table_05),colnames(gender))
anova_gender_gene_list_table_05<-merge(anova_gender_gene_list_table_05, gender, by=samecols05, all=TRUE)[samecols05]
a<-colnames(anova_gender_gene_list_table_05)
head(a)
anova_gender_gene_list_table_05<- anova_gender_gene_list_table_05[-c(1),]
nrow(gender)
nrow(anova_gender_gene_list_table_05)
head(colnames(anova_gender_gene_list_table_05), 30)
write.csv(anova_gender_gene_list_table_05, file = "/home/hwanghou/20200316/4838anova_gender_gene_list_table_05.csv", row.names = FALSE)



###  p<0.01

'''
gender<-tcga[,-c(1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)]
a<-colnames(gender)
head(a, 30)
table(gender$project_code_rename)
'''

#ttest_gender_bf<-read.csv("/home/solbi/20190118/tcga_4838_ttest_gender.csv")
head(ttest_gender_bf)
sapply(ttest_gender_bf, function(x) sum(is.na(x))) # count each columns 


anova_gender_p.B.01<-subset(ttest_gender_bf, P.Bonferroni<=0.01) # gene = 766
anova_gender_p.B.01 <- droplevels(anova_gender_p.B.01)
head(anova_gender_p.B.01)
str(anova_gender_p.B.01)

anova_gender_gene_list_01<-subset(anova_gender_p.B.01, select=c(rn))
colnames(anova_gender_gene_list_01) = c("gene_stable_id")
anova_gender_gene_list_01<-cbind(anova_gender_gene_list_01, donor_sex_rename=rep(1, nrow(anova_gender_gene_list_01)))
anova_gender_gene_list_01<-cbind(anova_gender_gene_list_01, num2=rep(2, nrow(anova_gender_gene_list_01)))
str(anova_gender_gene_list_01)
anova_gender_gene_list_01[, 2] <- as.numeric(as.character( anova_gender_gene_list_01[, 2] ))
anova_gender_gene_list_01[, 3] <- as.numeric(as.character( anova_gender_gene_list_01[, 3] ))
str(anova_gender_gene_list_01)

library(reshape2)
anova_gender_gene_list_table_01<-dcast(anova_gender_gene_list_01, donor_sex_rename ~ gene_stable_id, value.var = "num2", sum) #library(reshape2)
ncol(anova_gender_gene_list_table_01)
nrow(anova_gender_gene_list_table_01)
#write.csv(anova_gender_gene_list_table_01, file = "/home/solbi/20190118/4838anova_gender_gene_list_table_01.csv", row.names = FALSE)

#anova_gender_gene_list_table_01<-read.csv("/home/solbi/20190118/4838anova_gender_gene_list_table_01.csv")
str(anova_gender_gene_list_table_01)
samecols01 <- intersect(colnames(anova_gender_gene_list_table_01),colnames(gender))
anova_gender_gene_list_table_01<-merge(anova_gender_gene_list_table_01, gender, by=samecols01, all=TRUE)[samecols01]
a<-colnames(anova_gender_gene_list_table_01)
head(a)
anova_gender_gene_list_table_01<- anova_gender_gene_list_table_01[-c(1),]
nrow(gender)
nrow(anova_gender_gene_list_table_01)
head(colnames(anova_gender_gene_list_table_01), 30)
write.csv(anova_gender_gene_list_table_01, file = "/home/hwanghou/20200316/4838anova_gender_gene_list_table_01.csv", row.names = FALSE)


##############################################################################################################################################################
#########################    tcga all data anova about race    ##############################################################################################
#tcga<-read.csv("/home/solbi/20190118/tcga_primary_dataset_4838_geneexpres_dataset_4featuredata.csv")
a<-colnames(tcga)
head(a, 30)
tcga_all_race<-tcga[,-c(1, 2, 3, 4, 5, 6,7,8,9,10,12)]
a<-colnames(tcga_all_race)
head(a, 30)

library(dplyr)
ano_race<-t(sapply(tcga_all_race[-1], function(x) 
  +   unlist(anova(lm(x~race_rename, data=tcga_all_race)))))
library(data.table)
ano_race<-as.data.table(ano_race, keep.rownames=TRUE) #library(data.table)
ano_race_bf<-cbind(ano_race, P.Bonferroni = p.adjust(ano_race$`Pr(>F)1`, method="bonferroni"))
write.csv(ano_race_bf, file = "/home/hwanghou/20200316/tcga_4838_anova_race.csv", row.names = FALSE)
ano_race_bf<-read.csv("/home/solbi/20190118/tcga_4838_anova_race.csv")

sapply(ano_race_bf, function(x) sum(is.na(x))) # count each columns 
#rn          Df1          Df2      Sum.Sq1      Sum.Sq2     Mean.Sq1     Mean.Sq2     F.value1     F.value2      Pr..F.1 
#0            0            0            0            0            0            0           53        37889           53 
#Pr..F.2 P.Bonferroni 
#37889           53
summary(ano_race_bf$Pr..F.1)
summary(ano_race_bf$P.Bonferroni)
str(subset(ano_race_bf, P.Bonferroni==1)) # gene = 36446
ttest_race_p.B.05<-subset(ano_race_bf, P.Bonferroni<=0.05) # gene = 215
str(ttest_race_p.B.05)
ttest_race_p.B.05 <- droplevels(ttest_race_p.B.05)
head(ttest_race_p.B.05)

str(ttest_race_p.B.05)
ttest_race_p.B.01<-subset(ano_race_bf, P.Bonferroni<=0.01) # gene = 73
str(ttest_race_p.B.01)
ttest_race_p.B.005<-subset(ano_race_bf, P.Bonferroni<=0.005) # gene = 59
str(ttest_race_p.B.005)
ttest_race_p.B.001<-subset(ano_race_bf, P.Bonferroni<=0.001) # gene = 32
str(ttest_race_p.B.001)

str(subset(ano_race_bf, P.Bonferroni<=0.001))

###  p<0.05
race<-tcga[,-c(1, 2, 3, 4, 5, 6,7,8,9,10,12)]
a<-colnames(race)
head(a, 30)
table(race$race_rename)

ano_race_bf<-read.csv("/home/solbi/20190118/tcga_4838_anova_race.csv")
head(ano_race_bf)
sapply(ano_race_bf, function(x) sum(is.na(x))) # count each columns 


anova_race_p.B.05<-subset(ano_race_bf, P.Bonferroni<=0.05) # gene = 766
anova_race_p.B.05 <- droplevels(anova_race_p.B.05)
head(anova_race_p.B.05)
str(anova_race_p.B.05)

anova_race_gene_list_05<-subset(anova_race_p.B.05, select=c(rn))
colnames(anova_race_gene_list_05) = c("gene_stable_id")
anova_race_gene_list_05<-cbind(anova_race_gene_list_05, race_rename=rep(1, nrow(anova_race_gene_list_05)))
anova_race_gene_list_05<-cbind(anova_race_gene_list_05, num2=rep(2, nrow(anova_race_gene_list_05)))
str(anova_race_gene_list_05)
anova_race_gene_list_05[, 2] <- as.numeric(as.character( anova_race_gene_list_05[, 2] ))
anova_race_gene_list_05[, 3] <- as.numeric(as.character( anova_race_gene_list_05[, 3] ))
str(anova_race_gene_list_05)

library(reshape2)
anova_race_gene_list_table_05<-dcast(anova_race_gene_list_05, race_rename ~ gene_stable_id, value.var = "num2", sum) #library(reshape2)
ncol(anova_race_gene_list_table_05)
nrow(anova_race_gene_list_table_05)

str(anova_race_gene_list_table_05)
samecols05 <- intersect(colnames(anova_race_gene_list_table_05),colnames(race))
anova_race_gene_list_table_05<-merge(anova_race_gene_list_table_05, race, by=samecols05, all=TRUE)[samecols05]
a<-colnames(anova_race_gene_list_table_05)
head(a)
anova_race_gene_list_table_05<- anova_race_gene_list_table_05[-c(1),]
nrow(race)
nrow(anova_race_gene_list_table_05)
head(colnames(anova_race_gene_list_table_05), 10)
write.csv(anova_race_gene_list_table_05, file = "/home/hwanghou/20200316/4838anova_race_gene_list_table_05.csv", row.names = FALSE)


###  p<0.01
'''
race<-tcga[,-c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 14)]
a<-colnames(race)
head(a, 30)
table(race$project_code_rename)
'''

ano_race_bf<-read.csv("/home/solbi/20190118/tcga_4838_anova_race.csv")
head(ano_race_bf)
sapply(ano_race_bf, function(x) sum(is.na(x))) # count each columns 


anova_race_p.B.01<-subset(ano_race_bf, P.Bonferroni<=0.01) # gene = 766
anova_race_p.B.01 <- droplevels(anova_race_p.B.01)
head(anova_race_p.B.01)
str(anova_race_p.B.01)

anova_race_gene_list_01<-subset(anova_race_p.B.01, select=c(rn))
colnames(anova_race_gene_list_01) = c("gene_stable_id")
anova_race_gene_list_01<-cbind(anova_race_gene_list_01, race_rename=rep(1, nrow(anova_race_gene_list_01)))
anova_race_gene_list_01<-cbind(anova_race_gene_list_01, num2=rep(2, nrow(anova_race_gene_list_01)))
str(anova_race_gene_list_01)
anova_race_gene_list_01[, 2] <- as.numeric(as.character( anova_race_gene_list_01[, 2] ))
anova_race_gene_list_01[, 3] <- as.numeric(as.character( anova_race_gene_list_01[, 3] ))
str(anova_race_gene_list_01)

library(reshape2)
anova_race_gene_list_table_01<-dcast(anova_race_gene_list_01, race_rename ~ gene_stable_id, value.var = "num2", sum) #library(reshape2)
ncol(anova_race_gene_list_table_01)
nrow(anova_race_gene_list_table_01)
str(anova_race_gene_list_table_01)
samecols01 <- intersect(colnames(anova_race_gene_list_table_01),colnames(race))
anova_race_gene_list_table_01<-merge(anova_race_gene_list_table_01, race, by=samecols01, all=TRUE)[samecols01]
a<-colnames(anova_race_gene_list_table_01)
head(a)
anova_race_gene_list_table_01<- anova_race_gene_list_table_01[-c(1),]
nrow(race)
nrow(anova_race_gene_list_table_01)
head(colnames(anova_race_gene_list_table_01), 10)
write.csv(anova_race_gene_list_table_01, file = "/home/hwanghou/20200316/4838anova_race_gene_list_table_01.csv", row.names = FALSE)



##############################################################################################################################################################
#########################    tcga all data anova about age    ##############################################################################################
tcga<-read.csv("/home/solbi/20190118/tcga_primary_dataset_4838_geneexpres_dataset_4featuredata.csv")

tcga_all_age<-tcga[,-c(1, 2, 3, 4, 5, 6, 7, 9, 10, 11,12)]

a<-colnames(tcga_all_age)
head(a, 30)
tcga_all_age[-1]

library(dplyr)
ano_age<-t(sapply(tcga_all_age[-1], function(x) 
  +   unlist(anova(lm(x~donor_age_at_diagnosis_rename, data=tcga_all_age)))))
library(data.table)
ano_age<-as.data.table(ano_age, keep.rownames=TRUE) #library(data.table)
ano_age_bf<-cbind(ano_age, P.Bonferroni = p.adjust(ano_age$`Pr(>F)1`, method="bonferroni"))
write.csv(ano_age_bf, file = "/home/solbi/20190118/tcga_4838_anova_age.csv", row.names = FALSE)
ano_age_bf<-read.csv("/home/solbi/20190118/tcga_4838_anova_age.csv")

summary(ano_age_bf$Pr..F.1)
summary(ano_age_bf$P.Bonferroni)
str(subset(ano_age_bf, P.Bonferroni==1)) # gene = 36446
ttest_age_p.B.05<-subset(ano_age_bf, P.Bonferroni<=0.05) # gene = 30
str(ttest_age_p.B.05)
ttest_age_p.B.05 <- droplevels(ttest_age_p.B.05)
head(ttest_age_p.B.05)
str(ttest_age_p.B.05)
ttest_age_p.B.01<-subset(ano_age_bf, P.Bonferroni<=0.01) # gene = 24
str(ttest_age_p.B.01)
ttest_age_p.B.005<-subset(ano_age_bf, P.Bonferroni<=0.005) # gene = 20
str(ttest_age_p.B.005)
ttest_age_p.B.001<-subset(ano_age_bf, P.Bonferroni<=0.001) # gene = 19
str(ttest_age_p.B.001)

###  p<0.05
age<-tcga[,-c(1, 2, 3, 4, 5, 6, 7, 9, 10, 11,12)]
a<-colnames(age)
head(a, 30)
table(age$donor_age_at_diagnosis_rename)

ano_age_bf<-read.csv("/home/solbi/20190118/tcga_4838_anova_age.csv")
head(ano_age_bf)
sapply(ano_age_bf, function(x) sum(is.na(x))) # count each columns 


anova_age_p.B.05<-subset(ano_age_bf, P.Bonferroni<=0.05) # gene = 766
anova_age_p.B.05 <- droplevels(anova_age_p.B.05)
head(anova_age_p.B.05)
str(anova_age_p.B.05)

anova_age_gene_list_05<-subset(anova_age_p.B.05, select=c(rn))
colnames(anova_age_gene_list_05) = c("gene_stable_id")
anova_age_gene_list_05<-cbind(anova_age_gene_list_05, donor_age_at_diagnosis_rename=rep(1, nrow(anova_age_gene_list_05)))
anova_age_gene_list_05<-cbind(anova_age_gene_list_05, num2=rep(2, nrow(anova_age_gene_list_05)))
str(anova_age_gene_list_05)
anova_age_gene_list_05[, 2] <- as.numeric(as.character( anova_age_gene_list_05[, 2] ))
anova_age_gene_list_05[, 3] <- as.numeric(as.character( anova_age_gene_list_05[, 3] ))
str(anova_age_gene_list_05)

library(reshape2)
anova_age_gene_list_table_05<-dcast(anova_age_gene_list_05, donor_age_at_diagnosis_rename ~ gene_stable_id, value.var = "num2", sum) #library(reshape2)
ncol(anova_age_gene_list_table_05)
nrow(anova_age_gene_list_table_05)
str(anova_age_gene_list_table_05)
samecols05 <- intersect(colnames(anova_age_gene_list_table_05),colnames(age))
anova_age_gene_list_table_05<-merge(anova_age_gene_list_table_05, age, by=samecols05, all=TRUE)[samecols05]
a<-colnames(anova_age_gene_list_table_05)
head(a)
anova_age_gene_list_table_05<- anova_age_gene_list_table_05[-c(1),]
nrow(age)
nrow(anova_age_gene_list_table_05_2)
head(colnames(anova_age_gene_list_table_05), 10)
write.csv(anova_age_gene_list_table_05, file = "/home/hwanghou/20200316/4838anova_age_gene_list_table_05.csv", row.names = FALSE)


###  p<0.01
'''
age<-tcga[,-c(1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14)]
a<-colnames(age)
head(a, 30)
table(age$project_code_rename)
'''

ano_age_bf<-read.csv("/home/solbi/20190118/tcga_4838_anova_age.csv")
head(ano_age_bf)
sapply(ano_age_bf, function(x) sum(is.na(x))) # count each columns 


anova_age_p.B.01<-subset(ano_age_bf, P.Bonferroni<=0.01) # gene = 766
anova_age_p.B.01 <- droplevels(anova_age_p.B.01)
head(anova_age_p.B.01)
str(anova_age_p.B.01)

anova_age_gene_list_01<-subset(anova_age_p.B.01, select=c(rn))
colnames(anova_age_gene_list_01) = c("gene_stable_id")
anova_age_gene_list_01<-cbind(anova_age_gene_list_01, donor_age_at_diagnosis_rename=rep(1, nrow(anova_age_gene_list_01)))
anova_age_gene_list_01<-cbind(anova_age_gene_list_01, num2=rep(2, nrow(anova_age_gene_list_01)))
str(anova_age_gene_list_01)
anova_age_gene_list_01[, 2] <- as.numeric(as.character( anova_age_gene_list_01[, 2] ))
anova_age_gene_list_01[, 3] <- as.numeric(as.character( anova_age_gene_list_01[, 3] ))
str(anova_age_gene_list_01)

library(reshape2)
anova_age_gene_list_table_01<-dcast(anova_age_gene_list_01, donor_age_at_diagnosis_rename ~ gene_stable_id, value.var = "num2", sum) #library(reshape2)
ncol(anova_age_gene_list_table_01)
nrow(anova_age_gene_list_table_01)
write.csv(anova_age_gene_list_table_01, file = "/home/solbi/20190118/4838anova_age_gene_list_table_01.csv", row.names = FALSE)

anova_age_gene_list_table_01<-read.csv("/home/solbi/20190118/4838anova_age_gene_list_table_01.csv")
str(anova_age_gene_list_table_01)
samecols01 <- intersect(colnames(anova_age_gene_list_table_01),colnames(age))
anova_age_gene_list_table_01<-merge(anova_age_gene_list_table_01, age, by=samecols01, all=TRUE)[samecols01]
a<-colnames(anova_age_gene_list_table_01)
head(a)
anova_age_gene_list_table_01<- anova_age_gene_list_table_01[-c(1),]
nrow(age)
nrow(anova_age_gene_list_table_01)
head(colnames(anova_age_gene_list_table_01), 10)
write.csv(anova_age_gene_list_table_01, file = "/home/hwanghou/20200316/4838anova_age_gene_list_table_01.csv", row.names = FALSE)



##############################################################################################################################################################
#########################    tcga all data anova about cancer    ##############################################################################################
#tcga<-read.csv("/home/solbi/20190118/python_tcga_4feature.samples_rna.expression.data.csv")
tcga_all_cancer<-tcga[,-c(1, 2, 3, 4, 5, 6, 7, 8, 10,11,12)]
head(colnames(tcga),20)
a<-colnames(tcga_all_cancer)
head(a, 30)

library(dplyr)
ano_cancer<-t(sapply(tcga_all_cancer[-1], function(x) 
  +   unlist(anova(lm(x~project_code_rename, data=tcga_all_cancer)))))
library(data.table)
ano_cancer<-as.data.table(ano_cancer, keep.rownames=TRUE) #library(data.table)
ano_cancer_bf<-cbind(ano_cancer, P.Bonferroni = p.adjust(ano_cancer$`Pr(>F)1`, method="bonferroni"))
write.csv(ano_cancer_bf, file = "/home/solbi/20190118/tcga_4838_anova_cancer.csv", row.names = FALSE)
ano_cancer_bf<-read.csv("/home/solbi/20190118/tcga_4838_anova_cancer.csv")

summary(ano_cancer_bf$Pr..F.1)
summary(ano_cancer_bf$P.Bonferroni)
str(subset(ano_cancer_bf, P.Bonferroni==1)) # gene = 238
ttest_cancer_p.B.05<-subset(ano_cancer_bf, P.Bonferroni<=0.05) # gene = 37576
str(ttest_cancer_p.B.05)
ttest_cancer_p.B.05 <- droplevels(ttest_cancer_p.B.05)
head(ttest_cancer_p.B.05)
str(ttest_cancer_p.B.05)
ttest_cancer_p.B.01<-subset(ano_cancer_bf, P.Bonferroni<=0.01) # gene = 37572
str(ttest_cancer_p.B.01)
ttest_cancer_p.B.005<-subset(ano_cancer_bf, P.Bonferroni<=0.005) # gene = 37566
str(ttest_cancer_p.B.005)
ttest_cancer_p.B.001<-subset(ano_cancer_bf, P.Bonferroni<=0.001) # gene = 37560
str(ttest_cancer_p.B.001)

###  p<0.05
cancer<-tcga[,-c(1, 2, 3, 4, 5, 6, 7, 8, 10,11,12)]

a<-colnames(cancer)
head(a, 30)
table(cancer$project_code_rename)

ano_cancer_bf<-read.csv("/home/solbi/20190118/tcga_4838_anova_cancer.csv")
head(ano_cancer_bf)
sapply(ano_cancer_bf, function(x) sum(is.na(x))) # count each columns 


anova_cancer_p.B.05<-subset(ano_cancer_bf, P.Bonferroni<=0.05) # gene = 766
anova_cancer_p.B.05 <- droplevels(anova_cancer_p.B.05)
head(anova_cancer_p.B.05)
str(anova_cancer_p.B.05)

anova_cancer_gene_list_05<-subset(anova_cancer_p.B.05, select=c(rn))
colnames(anova_cancer_gene_list_05) = c("gene_stable_id")
anova_cancer_gene_list_05<-cbind(anova_cancer_gene_list_05, project_code_rename=rep(1, nrow(anova_cancer_gene_list_05)))
anova_cancer_gene_list_05<-cbind(anova_cancer_gene_list_05, num2=rep(2, nrow(anova_cancer_gene_list_05)))
str(anova_cancer_gene_list_05)
anova_cancer_gene_list_05[, 2] <- as.numeric(as.character( anova_cancer_gene_list_05[, 2] ))
anova_cancer_gene_list_05[, 3] <- as.numeric(as.character( anova_cancer_gene_list_05[, 3] ))
str(anova_cancer_gene_list_05)

library(reshape2)
anova_cancer_gene_list_table_05<-dcast(anova_cancer_gene_list_05, project_code_rename ~ gene_stable_id, value.var = "num2", sum) #library(reshape2)
ncol(anova_cancer_gene_list_table_05)
nrow(anova_cancer_gene_list_table_05)
write.csv(anova_cancer_gene_list_table_05, file = "/home/solbi/20190118/4838anova_cancer_gene_list_table_05.csv", row.names = FALSE)

anova_cancer_gene_list_table_05<-read.csv("/home/solbi/20190118/4838anova_cancer_gene_list_table_05.csv")
str(anova_cancer_gene_list_table_05)

samecols05 <- intersect(colnames(anova_cancer_gene_list_table_05),colnames(cancer))
anova_cancer_gene_list_table_05<-merge(anova_cancer_gene_list_table_05, cancer, by=samecols05, all=TRUE)[samecols05]
head(colnames(anova_cancer_gene_list_table_05))
anova_cancer_gene_list_table_05<- anova_cancer_gene_list_table_05[-c(1),]
nrow(cancer)
nrow(anova_cancer_gene_list_table_05)
write.csv(anova_cancer_gene_list_table_05, file = "/home/hwanghou/20200316/4838anova_cancer_gene_list_table_05.csv", row.names = FALSE)


###  p<0.01
'''
cancer<-tcga[,-c(1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14)]
a<-colnames(cancer)
head(a, 30)
table(cancer$project_code_rename)
'''

ano_cancer_bf<-read.csv("/home/solbi/20190118/tcga_4838_anova_cancer.csv")
head(ano_cancer_bf)
sapply(ano_cancer_bf, function(x) sum(is.na(x))) # count each columns 


anova_cancer_p.B.01<-subset(ano_cancer_bf, P.Bonferroni<=0.01) # gene = 766
anova_cancer_p.B.01 <- droplevels(anova_cancer_p.B.01)
head(anova_cancer_p.B.01)
str(anova_cancer_p.B.01)

anova_cancer_gene_list_01<-subset(anova_cancer_p.B.01, select=c(rn))
colnames(anova_cancer_gene_list_01) = c("gene_stable_id")
anova_cancer_gene_list_01<-cbind(anova_cancer_gene_list_01, project_code_rename=rep(1, nrow(anova_cancer_gene_list_01)))
anova_cancer_gene_list_01<-cbind(anova_cancer_gene_list_01, num2=rep(2, nrow(anova_cancer_gene_list_01)))
str(anova_cancer_gene_list_01)
anova_cancer_gene_list_01[, 2] <- as.numeric(as.character( anova_cancer_gene_list_01[, 2] ))
anova_cancer_gene_list_01[, 3] <- as.numeric(as.character( anova_cancer_gene_list_01[, 3] ))
str(anova_cancer_gene_list_01)

library(reshape2)
anova_cancer_gene_list_table_01<-dcast(anova_cancer_gene_list_01, project_code_rename ~ gene_stable_id, value.var = "num2", sum) #library(reshape2)
ncol(anova_cancer_gene_list_table_01)
nrow(anova_cancer_gene_list_table_01)
#write.csv(anova_cancer_gene_list_table_01, file = "/home/solbi/20190118/4838anova_cancer_gene_list_table_01.csv", row.names = FALSE)
#anova_cancer_gene_list_table_01<-read.csv("/home/solbi/20190118/4838anova_cancer_gene_list_table_01.csv")
str(anova_cancer_gene_list_table_01)
samecols01 <- intersect(colnames(anova_cancer_gene_list_table_01),colnames(cancer))
anova_cancer_gene_list_table_01<-merge(anova_cancer_gene_list_table_01, cancer, by=samecols01, all=TRUE)[samecols01]
a<-colnames(anova_cancer_gene_list_table_01)
head(a)
anova_cancer_gene_list_table_01<- anova_cancer_gene_list_table_01[-c(1),]
nrow(cancer)
nrow(anova_cancer_gene_list_table_01)
head(colnames(anova_cancer_gene_list_table_01), 10)
write.csv(anova_cancer_gene_list_table_01, file = "/home/hwanghou/20200316/4838anova_cancer_gene_list_table_01.csv", row.names = FALSE)


##############################################################################################################################################################
#########################    tcga all data anova about stage    ##############################################################################################
tcga5477<-read.csv("/home/solbi/20190118/tcga_primary_dataset_3782_geneexpres_dataset_stage_featuredata.csv")
tcga_all_stage<-tcga_primary_dataset_4838_geneexpres_dataset5[,-c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11)]
a<-colnames(tcga_all_stage)
head(a, 30)

library(dplyr)
ano_stage<-t(sapply(tcga_all_stage[-1], function(x) 
  +   unlist(anova(lm(x~Stage_rename, data=tcga_all_stage)))))
library(data.table)
ano_stage<-as.data.table(ano_stage, keep.rownames=TRUE) #library(data.table)
ano_stage_bf<-cbind(ano_stage, P.Bonferroni = p.adjust(ano_stage$`Pr(>F)1`, method="bonferroni"))
write.csv(ano_stage_bf, file = "/home/solbi/20190118/tcga_3782_anova_stage.csv", row.names = FALSE)
ano_stage_bf<-read.csv("/home/solbi/20190118/tcga_3782_anova_stage.csv")

summary(ano_stage_bf$Pr..F.1)
summary(ano_stage_bf$P.Bonferroni)
str(subset(ano_stage_bf, P.Bonferroni==1)) 
ttest_stage_p.B.05<-subset(ano_stage_bf, P.Bonferroni<=0.05)
str(ttest_stage_p.B.05)
ttest_stage_p.B.05 <- droplevels(ttest_stage_p.B.05)
head(ttest_stage_p.B.05)
str(ttest_stage_p.B.05)
str(subset(ano_stage_bf, P.Bonferroni<=0.05)) #3129
str(subset(ano_stage_bf, P.Bonferroni<=0.01)) # gene = 2673
str(subset(ano_stage_bf, P.Bonferroni<=0.005)) # gene = 2465
str(subset(ano_stage_bf, P.Bonferroni<=0.001)) # gene = 2113

###  p<0.05
stage<-tcga_primary_dataset_4838_geneexpres_dataset5[,-c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11)]
head(colnames(tcga_primary_dataset_4838_geneexpres_dataset5),20)
a<-colnames(stage)
head(a, 30)
table(stage$Stage_rename)

ano_stage_bf<-read.csv("/home/solbi/20190118/tcga_5477_anova_stage.csv")
head(ano_stage_bf)
sapply(ano_stage_bf, function(x) sum(is.na(x))) # count each columns 


anova_stage_p.B.05<-subset(ano_stage_bf, P.Bonferroni<=0.05) # gene = 766
anova_stage_p.B.05 <- droplevels(anova_stage_p.B.05)
head(anova_stage_p.B.05)
str(anova_stage_p.B.05)

anova_stage_gene_list_05<-subset(anova_stage_p.B.05, select=c(rn))
colnames(anova_stage_gene_list_05) = c("gene_stable_id")
anova_stage_gene_list_05<-cbind(anova_stage_gene_list_05, Stage_rename=rep(1, nrow(anova_stage_gene_list_05)))
anova_stage_gene_list_05<-cbind(anova_stage_gene_list_05, num2=rep(2, nrow(anova_stage_gene_list_05)))
str(anova_stage_gene_list_05)
anova_stage_gene_list_05[, 2] <- as.numeric(as.character( anova_stage_gene_list_05[, 2] ))
anova_stage_gene_list_05[, 3] <- as.numeric(as.character( anova_stage_gene_list_05[, 3] ))
str(anova_stage_gene_list_05)

library(reshape2)
anova_stage_gene_list_table_05<-dcast(anova_stage_gene_list_05, Stage_rename ~ gene_stable_id, value.var = "num2", sum) #library(reshape2)
ncol(anova_stage_gene_list_table_05)
nrow(anova_stage_gene_list_table_05)
str(anova_stage_gene_list_table_05)
samecols05 <- intersect(colnames(anova_stage_gene_list_table_05),colnames(stage))
anova_stage_gene_list_table_05<-merge(anova_stage_gene_list_table_05, stage, by=samecols05, all=TRUE)[samecols05]
head(colnames(anova_stage_gene_list_table_05))
anova_stage_gene_list_table_05<- anova_stage_gene_list_table_05[-c(1),]
nrow(stage)
nrow(anova_stage_gene_list_table_05)
head(colnames(anova_stage_gene_list_table_05), 10)
write.csv(anova_stage_gene_list_table_05, file = "/home/hwanghou/20200316/3782anova_stage_gene_list_table_05.csv", row.names = FALSE)


###  p<0.01
'''
stage<-tcga[,-c(1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)]
a<-colnames(stage)
head(a, 30)
table(stage$project_code_rename)
'''

ano_stage_bf<-read.csv("/home/solbi/20190118/tcga_5477_anova_stage.csv")
head(ano_stage_bf)
sapply(ano_stage_bf, function(x) sum(is.na(x))) # count each columns 


anova_stage_p.B.01<-subset(ano_stage_bf, P.Bonferroni<=0.01) # gene = 766
anova_stage_p.B.01 <- droplevels(anova_stage_p.B.01)
head(anova_stage_p.B.01)
str(anova_stage_p.B.01)

anova_stage_gene_list_01<-subset(anova_stage_p.B.01, select=c(rn))
colnames(anova_stage_gene_list_01) = c("gene_stable_id")
anova_stage_gene_list_01<-cbind(anova_stage_gene_list_01, Stage_rename=rep(1, nrow(anova_stage_gene_list_01)))
anova_stage_gene_list_01<-cbind(anova_stage_gene_list_01, num2=rep(2, nrow(anova_stage_gene_list_01)))
str(anova_stage_gene_list_01)
anova_stage_gene_list_01[, 2] <- as.numeric(as.character( anova_stage_gene_list_01[, 2] ))
anova_stage_gene_list_01[, 3] <- as.numeric(as.character( anova_stage_gene_list_01[, 3] ))
str(anova_stage_gene_list_01)

library(reshape2)
anova_stage_gene_list_table_01<-dcast(anova_stage_gene_list_01, Stage_rename ~ gene_stable_id, value.var = "num2", sum) #library(reshape2)
ncol(anova_stage_gene_list_table_01)
nrow(anova_stage_gene_list_table_01)
str(anova_stage_gene_list_table_01)
samecols01 <- intersect(colnames(anova_stage_gene_list_table_01),colnames(stage))
anova_stage_gene_list_table_01<-merge(anova_stage_gene_list_table_01, stage, by=samecols01, all=TRUE)[samecols01]
a<-colnames(anova_stage_gene_list_table_01)
head(a)
anova_stage_gene_list_table_01<- anova_stage_gene_list_table_01[-c(1),]
nrow(stage)
nrow(anova_stage_gene_list_table_01)
head(colnames(anova_stage_gene_list_table_01), 10)
head(colnames(anova_stage_gene_list_table_01), 10)
write.csv(anova_stage_gene_list_table_01, file = "/home/hwanghou/20200316/3782anova_stage_gene_list_table_01.csv", row.names = FALSE)

###################################################################################################
###################################################################################################
#####################################   top 5 cacner type         #################################
###################################################################################################
###################################################################################################
### save 
tcga<-read.csv("/home/solbi/20190118/tcga_primary_dataset_4838_geneexpres_dataset_4featuredata.csv")
a<-split(tcga, tcga$project_code)
lapply(1:length(a), function (x) write.csv(a[[x]],file = paste('/home/solbi/20190118/cancersample/tcga_primary_dataset_4838_geneexpres_dataset_4featuredata_', names (a[x]),'.csv',sep=""), row.names = F))
### 
### save 
tcga3782<-read.csv("/home/solbi/20190118/tcga_primary_dataset_3782_geneexpres_dataset_stage_featuredata.csv")
a<-split(tcga3782, tcga3782$project_code)
lapply(1:length(a), function (x) write.csv(a[[x]],file = paste('/home/solbi/20190118/cancersample/tcga_primary_dataset_3782_geneexpres_dataset_stage_featuredata_', names (a[x]),'.csv',sep=""), row.names = F))


################################################################################################################################################################################################
############################  1. BRCA data  ######################################################################################################################################################
BRCA<-read.csv("/home/hwanghou/20200316/cancersample/tcga_primary_dataset_3782_geneexpres_dataset_stage_featuredata_BRCA-US.csv")

table(BRCA$donor_sex_rename)

table(BRCA$donor_age_at_diagnosis_rename)
#1   2   3   4   5   6   7   8 
#7  58 176 229 213 113  44   7 
BRCA_python<-BRCA
BRCA_python$donor_age_at_diagnosis_rename<-gsub("1", "0", BRCA_python$donor_age_at_diagnosis_rename)
BRCA_python$donor_age_at_diagnosis_rename<-gsub("2", "1", BRCA_python$donor_age_at_diagnosis_rename)
BRCA_python$donor_age_at_diagnosis_rename<-gsub("3", "2", BRCA_python$donor_age_at_diagnosis_rename)
BRCA_python$donor_age_at_diagnosis_rename<-gsub("4", "3", BRCA_python$donor_age_at_diagnosis_rename)
BRCA_python$donor_age_at_diagnosis_rename<-gsub("5", "4", BRCA_python$donor_age_at_diagnosis_rename)
BRCA_python$donor_age_at_diagnosis_rename<-gsub("6", "5", BRCA_python$donor_age_at_diagnosis_rename)
BRCA_python$donor_age_at_diagnosis_rename<-gsub("7", "6", BRCA_python$donor_age_at_diagnosis_rename)
BRCA_python$donor_age_at_diagnosis_rename<-gsub("8", "7", BRCA_python$donor_age_at_diagnosis_rename)
table(BRCA_python$donor_age_at_diagnosis_rename)

table(BRCA_python$race_rename)
BRCA_python$race_rename<-gsub("4", "3", BRCA_python$race_rename)
table(BRCA_python$race)
head(colnames(BRCA_python), 30)

write.csv(BRCA_python, file = "/home/hwanghou/20200316/cancersample/tcga_primary_dataset_3782_geneexpres_dataset_stage_featuredata_BRCA-US_python.csv", row.names = FALSE)

BRCA_python<-read.csv("/home/solbi/20190118/cancersample/tcga_primary_dataset_4838_geneexpres_dataset_4featuredata_BRCA-US.csv")


################################################################################################################################################################################################
############################  2. KIRC data  ######################################################################################################################################################
KIRC<-read.csv("/home/hwanghou/20200316/cancersample/tcga_primary_dataset_3782_geneexpres_dataset_stage_featuredata_KIRC-US.csv")

table(KIRC$donor_sex_rename)

table(KIRC$donor_age_at_diagnosis_rename)
#1   2   3   4   5   6   7   8 
#7  58 176 229 213 113  44   7 
KIRC_python<-KIRC
KIRC_python$donor_age_at_diagnosis_rename<-gsub("1", "0", KIRC_python$donor_age_at_diagnosis_rename)
KIRC_python$donor_age_at_diagnosis_rename<-gsub("2", "1", KIRC_python$donor_age_at_diagnosis_rename)
KIRC_python$donor_age_at_diagnosis_rename<-gsub("3", "2", KIRC_python$donor_age_at_diagnosis_rename)
KIRC_python$donor_age_at_diagnosis_rename<-gsub("4", "3", KIRC_python$donor_age_at_diagnosis_rename)
KIRC_python$donor_age_at_diagnosis_rename<-gsub("5", "4", KIRC_python$donor_age_at_diagnosis_rename)
KIRC_python$donor_age_at_diagnosis_rename<-gsub("6", "5", KIRC_python$donor_age_at_diagnosis_rename)
KIRC_python$donor_age_at_diagnosis_rename<-gsub("7", "6", KIRC_python$donor_age_at_diagnosis_rename)
KIRC_python$donor_age_at_diagnosis_rename<-gsub("8", "7", KIRC_python$donor_age_at_diagnosis_rename)
table(KIRC_python$donor_age_at_diagnosis_rename)
table(KIRC_python$Stage_rename)

table(KIRC_python$race_rename)
KIRC_python$race_rename<-gsub("1", "0", KIRC_python$race_rename)
KIRC_python$race_rename<-gsub("2", "1", KIRC_python$race_rename)
KIRC_python$race_rename<-gsub("4", "2", KIRC_python$race_rename)
table(KIRC_python$race)
head(colnames(KIRC_python), 30)

write.csv(KIRC_python, file = "/home/hwanghou/20200316/cancersample/tcga_primary_dataset_3782_geneexpres_dataset_stage_featuredata_KIRC-US_python.csv", row.names = FALSE)


################################################################################################################################################################################################
############################  3. UCEC data  ######################################################################################################################################################
UCEC<-read.csv("/home/hwanghou/20200316/cancersample/tcga_primary_dataset_4838_geneexpres_dataset_4featuredata_UCEC-US.csv")

table(UCEC$donor_sex_rename)

table(UCEC$donor_age_at_diagnosis_rename)
# 2   3   4   5   6   7   8
#13  25 111 164  92  35   3

UCEC_python<-UCEC
UCEC_python$donor_age_at_diagnosis_rename<-gsub("2", "0", UCEC_python$donor_age_at_diagnosis_rename)
UCEC_python$donor_age_at_diagnosis_rename<-gsub("3", "1", UCEC_python$donor_age_at_diagnosis_rename)
UCEC_python$donor_age_at_diagnosis_rename<-gsub("4", "2", UCEC_python$donor_age_at_diagnosis_rename)
UCEC_python$donor_age_at_diagnosis_rename<-gsub("5", "3", UCEC_python$donor_age_at_diagnosis_rename)
UCEC_python$donor_age_at_diagnosis_rename<-gsub("6", "4", UCEC_python$donor_age_at_diagnosis_rename)
UCEC_python$donor_age_at_diagnosis_rename<-gsub("7", "5", UCEC_python$donor_age_at_diagnosis_rename)
UCEC_python$donor_age_at_diagnosis_rename<-gsub("8", "6", UCEC_python$donor_age_at_diagnosis_rename)
table(UCEC_python$donor_age_at_diagnosis_rename)

table(UCEC_python$race_rename)
#0   1   2   3   4
#4  19  63   8 349

write.csv(UCEC_python, file = "/home/hwanghou/20200316/cancersample/tcga_primary_dataset_4838_geneexpres_dataset_4featuredata_UCEC-US_python.csv", row.names = FALSE)

################################################################################################################################################################################################
############################  4. LUAD data  ######################################################################################################################################################
LUAD<-read.csv("/home/hwanghou/20200316/cancersample/tcga_primary_dataset_3782_geneexpres_dataset_stage_featuredata_LUAD-US.csv")

table(LUAD$donor_sex_rename)

table(LUAD$donor_age_at_diagnosis_rename)
#  2   3   4   5   6   7 
#  1  20  81 121 127  26 
LUAD_python<-LUAD
LUAD_python$donor_age_at_diagnosis_rename<-gsub("2", "0", LUAD_python$donor_age_at_diagnosis_rename)
LUAD_python$donor_age_at_diagnosis_rename<-gsub("3", "1", LUAD_python$donor_age_at_diagnosis_rename)
LUAD_python$donor_age_at_diagnosis_rename<-gsub("4", "2", LUAD_python$donor_age_at_diagnosis_rename)
LUAD_python$donor_age_at_diagnosis_rename<-gsub("5", "3", LUAD_python$donor_age_at_diagnosis_rename)
LUAD_python$donor_age_at_diagnosis_rename<-gsub("6", "4", LUAD_python$donor_age_at_diagnosis_rename)
LUAD_python$donor_age_at_diagnosis_rename<-gsub("7", "5", LUAD_python$donor_age_at_diagnosis_rename)
table(LUAD_python$donor_age_at_diagnosis_rename)

table(LUAD_python$race_rename)
LUAD_python$race_rename<-gsub("4", "3", LUAD_python$race_rename)
table(LUAD_python$race_rename)

write.csv(LUAD_python, file = "/home/hwanghou/20200316/cancersample/tcga_primary_dataset_3782_geneexpres_dataset_stage_featuredata_LUAD-US_python.csv", row.names = FALSE)



################################################################################################################################################################################################
############################  5. THCA data  ######################################################################################################################################################
THCA<-read.csv("/home/hwanghou/20200316/cancersample/tcga_primary_dataset_3782_geneexpres_dataset_stage_featuredata_THCA-US.csv")

table(THCA$donor_sex_rename)

table(THCA$donor_age_at_diagnosis_rename)

table(THCA$race_rename)
THCA_python<-THCA
THCA_python$race_rename<-gsub("4", "3", THCA_python$race_rename)
table(THCA_python$race_rename)

write.csv(THCA_python, file = "/home/hwanghou/20200316/cancersample/tcga_primary_dataset_3782_geneexpres_dataset_stage_featuredata_THCA-US_python.csv", row.names = FALSE)



