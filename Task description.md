### Пайплайн
Общая логика работы пайплайна такова:

1. Референсный геном подготавливается для анализа
1. С использованием референсного генома и файлов с данными проводится поиск генетических вариантов.  Для поиска генетических вариантов используется три различных инструмента: freebayes, samtools и HaplotypeCaller.
1. Производится подсчет генетических вариантов в файлах, полученных каждым из инструментов.
1. Производится поиск вариантов, общих для файлов, полученных каждым из инструментов.
1. Создается файл с отчетом (см. ниже)




#### Входные данные
В директории input размещены файлы, соответствующие результатам секвенирования. Это BAM файлы с расширением bam, каждый файл соответствует одному образцу.
В директории data размещен файл 22.fa – это файл референса.
Выходные данные
Выходные данные – текстовые файлы, размещаемые в директории output/reports. Файл имеет следующую структуру (в качестве разделителя используется tab):

vc_1 vc_1&vc_2 vc_1&vc_3
vc_2 None vc2&vc_3
vc_3 None None

vc_N – количество вариантов, полученных при использовании variant caller'а N (целое число)
vc_N&vc_M – количество вариантов, общих для variant caller'ов N и M (целое число)
variant caller'ы нумеруются в алфавитном порядке: freebayes, haplotypecaller, samtools
Операции и команды
В процессе работы вам потребуются следующие инструменты и команды:

Reference index
Создание индекса для файла референса, используется samtools:

samtools faidx data/22.fa 
В результате формируется файл data/22.fa.fai, которые требуется для работы variant caller'ов.

Reference dict
Формирование словаря контигов, используется picard. В тестовом и проверяжщем контейнере picard.jar доступен с использованием переменной окружения PICARD:

java -jar $PICARD CreateSequenceDictionary R=data/22.fa O=data/22.dict
Index BAM
Индексирования файла с исходными данными требуется для работы многих инструментов, в том числе и variant caller'ов.

samtools index input/sample.bam

В результате формируется файл input/sample.bam.bai
Variant Calling HaplotypeCaller
Поиск генетических вариантов с использованием HaplotypeCaller'а. В тестовом и проверяющем контейнерах GenomeAnalysisTK.jar доступен с использованием переменной окружения GATK:

java -jar $GATK -R data/22.fa -T HaplotypeCaller -I input/sample.bam -o output/any-path-you-like/sample.vcf

В результате формируется файл VCF, путь и имя вы задаете самостоятельно.
Variant Calling Freebayes
Поиск генетических вариантов с использованием freebayes:

freebayes -f data/22.fa input/sample.bam > output/any-path-you-like/sample.vcf

В результате формируется файл VCF, путь и имя вы задаете самостоятельно. 
Variant Calling Samtools
Поиск генетических вариантов с использованием samtools:
samtools mpileup -uf data/22.fa input/sample.bam | bcftools view -vcg - > output/any-path-you-like/sample.vcf

В результате формируется файл VCF, путь и имя вы задаете самостоятельно. 
Index Variant Files
Для поиска вариантов, общих для двух файлов, файлы необходимо проиндексировать, используются bgzip и tabix:

bgzip sample.vcf
tabix -p vcf sample.vcf.gz

Intersect Variants
Отбор вариантов, общих для двух файлов, используется vcf-isec из пакета vcftools:

vcf-isec -f -n +2 file_1.vcf.gz file_2.vcf.gz | bgzip -c > variants_common_for_file_1_and_file_2.vcf.gz
Общие варианты попадают в результирующий VCF файл, имя можно задать самостоятельно.

Важно: разные инструменты генерируют файлы вариантов, которые отличаются друг от друга, поэтому при сравнении мы используем ключ -f (сравнить, несмотря на несовпадение версий). При таком сравнении порядок файлов имеет значение! Для получения правильного ответа нужно использовать при сравнении пары: freebayes / haplotypecaller,  freebayes / samtools,  haplotypecaller / samtools.

Count Variants
Для подсчета вариантов в VCF файлах используется vcftools:

vcftools --vcf output/sample.vcf

В stderr выводится сообщение, включающее в себя информацию о количестве вариантов в файле:
VCFtools - v0.1.12
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
  --vcf sample.vcf

Using zlib version: 1.2.8
After filtering, kept 1 out of 1 Individuals
After filtering, kept 1058 out of a possible 1058 Sites
Run Time = 0.00 seconds
 
В этом примере количество вариантов в файле – 1058.

Для подсчета вариантов в сжатых VCF файлах используется аналогичный подход:
vcftools --gzvcf sample.vcf.gz

Report Creation / All
В результате для каждого образца из папки input должен быть сформирован отчет в папке output/reports. Формат отчета описан выше.
