# **Identification of genes associated with immune response using open database ImmuneSpace**

**Study SDY984 (immune response after varicella zoster vaccine)**

author: Мироненко Ольга

## **Общая информация**

<br>

В данном исследовательском проекте будут представлены результаты анализа данных по экспрессии генов среди добровольцев, прошедших **вакцинацию от varicella zoster живой ослабленной вакциной Zostavax**. Данные были получены в рамках исследования [SDY984](https://www.immport.org/shared/study/SDY984) по программе [Human Immunology Project Consortium (HIPC)](https://www.immunespace.org/) и выложены на страничке [HIPC-II Immune Signatures Data Resource and Analysis (IS2)](https://datatools.immunespace.org/project/home/Integrative_Public_Study/begin.view?SDY=IS2) (доступны для скачивания после регистрации). Результаты оригинального исследования по полной базе данных были опубликованы в статье (Li et al. 2017).

**Задача нашего исследовательского проекта** - определить, динамика экспрессии каких именно генов связана с иммунным ответом на вакцинацию и в какие сигнальные пути входят продукты этих генов. Первую часть этой задачи можно решать двояко: 

- с помощью оценки различий в экспрессии генов между испытуемыми с разным уровнем ответа на вакцинацию, 

- с помощью оценки взаимосвязи между экспрессией генов (и её динамикой) и итоговым ответом на вакцинацию.

<br>

**Процесс генерации иммунного ответа на вакцинацию** может быть схематично описан следующим образом (Pollard and Bijker 2021):

<img src="../OlgaMironenko/images/ImmResp_Vacc1.jpg" width="100%" />

<br>

Или так (Desmet and Ishii 2012):

<img src="../OlgaMironenko/images/ImmResp_Vacc2.jpg" width="100%" />

<br>

Вакцина может содержать молекулярные паттерны, ассоциированные с патогеном (PAMPs), или спровоцировать локальную реакцию в виде высвобождения молекулярных паттернов, ассоциированных с повреждением (DAMPs). Эти паттерны распознаются рецепторами распознавания образов (PRRs), экспрессируемых клетками, являющимися частью системы **врождённого иммунитета**: макрофагами, тучными клетками, нейтрофилами, дендритными клетками (DCs). Эти клетки призваны не определить, с чем именно - с какой именно инфекцией или иным патогеном - столкнулся организм, а, в зависимости от того, о каких клетках идёт речь, уничтожить или поглотить патоген, расщепить его на фрагменты и/ или просигнализировать клеткам адаптивного иммунитета о том, что организм с ним столкнулся. Система **адаптивного иммунитета** уже будет работать на то, чтобы сформировать/ задействовать специфичный, прицельный иммунный ответ именно на данный патоген, в частности, в виде выработки специфических антител к нему. 

Также в рамках врождённого иммунитета PAMPs и DAMPs могут распознаваться и лимфоидными клетками (например, NK-клетками), которые в ответ на это продуцируют цитокины, участвующие, в том числе, в активации и ориентации дендритных клеток, играющих основную роль в передаче сигнала от системы врождённого к системе адаптивного иммунного ответа. После активации дендритные клетки мигрируют в лимфатические узлы, где презентируют фрагменты поглощенного антигена T-клеткам через специальный механизм костимуляции рецепторов T-клеток (TCR), в процессе которой одни T-клетки (T-киллеры, или CD8+ лимфоциты) взаимодействуют с молекулами главного комплекса гистосовместимости первого типа (MCH-I) и при обнаружении чужеродных MHC-I выделяют белки, растворяющие патоген, а другие (T-хелперы, или CD4+ лимфоциты) взаимодействуют с молекулами второго типа (MHC-II) и при обранужении чужеродных MCH-II продуцируют хемокины и цитокины. В зависимости от цитокиновой среды CD4+ лимфоциты могут дифференцироваться в различные подтипы T-хелперов (Th), а также T-фолликулярные хелперы (Tfh) - последние активируют B-клетки, в процессе дифференцировки которых образуются плазматические клетки, вырабатывающие антиген-специфические антитела. В процессе дифференцировки B-клеток и CD8+ лимфоцитов также будут образовываться клетки памяти, что позволит иммунной системе в дальнейшем быстрее распознавать и реагировать на контакт с уже знакомым патогеном.

Описанная схема формирования иммунного ответа даёт примерное представление о том, какие биологические процессы она задействует. 

<br>

Что касается участия продуктов определенных генов в этих процессах, то представление о нём можно получить из **[Gene Ontology (GO)](http://geneontology.org/)** - специальной базы знаний, в которой хранится наиболее актуальная и полная на сегодняшний день информация в данной области, а именно: аннотация о молекулярных функциях белковых продуктов генов, биологических процессах, частью которых являются эти функции, и компонентах клеток, в которых продукты гена эти функции выполняют. В нашем исследовании нас будут интересовать, прежде всего, биологические процессы, в которых участвуют продукты генов-находок (генов с дифференциальной экспрессией для групп испытуемых с разной силой ответа на вакцинацию или генов, экспрессия которых окажется значимым образом связанной с ответом).

Дополнительно мы будем выделять те **сигнальные пути, которые непосредственно связаны с формированием иммунного ответа** на вакцинацию, а именно: процессы, являющиеся потомками следующих биологических процессов (на рисунках ниже показаны только "дети" соответствующих процессов):

- immune system process ([GO:0002376](https://www.ebi.ac.uk/QuickGO/term/GO:0002376)) - any process involved in the development or functioning of the immune system, an organismal system for calibrated responses to potential internal or invasive threats.

<a href="https://raw.githubusercontent.com/ejeej/ImmuneSpace_SDY984/main/OlgaMironenko/images/GO_0002376.png">
<img src="../OlgaMironenko/images/GO_0002376.png"/>
</a>

- defense response ([GO:0006952](https://www.ebi.ac.uk/QuickGO/term/GO:0006952)) - reactions, triggered in response to the presence of a foreign body or the occurrence of an injury, which result in restriction of damage to the organism attacked or prevention/recovery from the infection caused by the attack.  

<a href="https://raw.githubusercontent.com/ejeej/ImmuneSpace_SDY984/main/OlgaMironenko/images/GO_0006952.png">
<img src="../OlgaMironenko/images/GO_0006952.png"/>
</a>

- cytokine production ([GO:0001816](https://www.ebi.ac.uk/QuickGO/term/GO:0001816)) - appearance of a cytokine due to biosynthesis or secretion following a cellular stimulus, resulting in an increase in its intracellular or extracellular levels. 

<a href="https://raw.githubusercontent.com/ejeej/ImmuneSpace_SDY984/main/OlgaMironenko/images/GO_0001816.png">
<img src="../OlgaMironenko/images/GO_0001816.png"/>
</a>

- response to cytokine ([GO:0034097](https://www.ebi.ac.uk/QuickGO/term/GO:0034097)) - any process that results in a change in state or activity of a cell or an organism (in terms of movement, secretion, enzyme production, gene expression, etc.) as a result of a cytokine stimulus.

<a href="https://raw.githubusercontent.com/ejeej/ImmuneSpace_SDY984/main/OlgaMironenko/images/GO_0034097.png">
<img src="../OlgaMironenko/images/GO_0034097.png"/>
</a>

- MHC protein complex assembly ([GO:0002396](https://www.ebi.ac.uk/QuickGO/term/GO:0002396)) - the aggregation, arrangement and bonding together of a set of components to form an MHC protein complex.

<center>
<a href="https://raw.githubusercontent.com/ejeej/ImmuneSpace_SDY984/main/OlgaMironenko/images/GO_0002396.png">
<img src="../OlgaMironenko/images/GO_0002396.png" style="float:center;width:30%"/>
</a>
</center>

<br>

Понятно, что иммунный ответ задействует и другие сигнальные пути, помимо обозначенных, например, процессы, связанные с активацией клеток, их диффференцировкой, межклеточными взаимодействиями и т.п., но эти процессы могут происходить в организме и не в связи с вакцинацией, поэтому мы будем обращать на них внимание, но не будем относить к группе непосредственно связанных с иммуным ответом.

<br>

В выбранном нами исследовании **экспрессия генов** оценивалась у 35 добровольцев до вакцинации и в нескольких точках (1, 3, 7 дней) после неё. 

Также по завершении исследования (через 30 дней от начала) оценивался **ответ на вакцинацию** с разбиением испытуемых на 3 группы: low, moderate и high responder. Разбиение проводилось отдельно по нескольким критериям, но мы для определения статуса по ответу на вакцинацию будем использовать только MFC_p40. MFC - это maximum fold change (максимальная кратность изменения) для титра антител после вакцинации по сравнению с состоянием до неё (максимальная - поскольку некоторые вакцины в исследовании содержали несколько штаммов вируса и/или титр антител для них определялся с помощью нескольких тестов (Fourati et al. 2022)). В исследовании, с данными которого мы будем работать, штамм был один и титр антител измерялся с помощью одной методики - оценки уровня IgG методом ELISA (enzyme-linked immunosorbent assay, или иммуноферментный анализ, ИФА), поэтому в данном случае MFC фактически равна просто кратности изменения титра антител (формально это $log_2(FC)$). Что касается разбивки на группы по ответу, то оно производилось по перцентилям полученного значения для MFC (Fourati et al. 2022), а именно: для MFC_p40 все испытуемые, у кого MFC был меньше или равен 40-му перцентилю, считались low responder'ами, с MFC, равным или большим 60-го перцентиля - high responder'ами, остальные - moderate responder'ами (детали можно найти, в том числе, в [функциях](https://rdrr.io/github/RGLab/ImmuneSignatures2/src/R/immuneResponseCallGeneration.R), которые использовались для получения того датасета, с которым мы работаем).

В целях упрощения в рамках нашего исследовательского проекта мы исключим из анализа испытуемых с moderate ответом на вакцинацию (7 чел.), а также, на начальном этапе, не будем делать разбивку в зависимости от ветви исследования (возрастной группы испытуемых).

<br>

В **оригинальном исследовании** (Li et al. 2017) использовался больший по сравнению с доступным нам объём данных (и по количеству испытуемых, и по количеству показателей, в том числе касающихся метаболомики и проточной цитометрии, и по числу точек оценки последних) и решался более широкий спектр задач с несколько отличающимся от нашего акцентом, а именно: оценка взаимосвязи между изменениями метаболических показателей с изменениями экспрессии генов и их отдельных модулей на протяжении периода до 180 дней после вакцинации. Соответственно, авторам удалось более детально описать механизм ответа на рассматриваемую вакцину, без его дихотомизации. В своём исследовательском проекте мы ставим более узкую задачу - выявить гены, дифференциально экспрессированные у испытуемых с сильным и слабым ответом на вакцинацию или экспрессия которых связана (положительно или отрицательно) с вероятностью сильного ответа, а также те сигнальные пути, в которых задйствованы продукты этих генов.

<br>

## **Описательная статистика: характеристики испытуемых, baseline**

<br>

Описательная статистика по участникам исследования на начало исследования (до вакцинации) представлена в Table 1 ниже. Поскольку мы планируем в своём исследовании, в том числе, сравнивать испытуемых с разным уровнем ответа по экспрессии генов, данные в таблице представлены в разбивке по этим группам. Young - добровольцы в возрасте 25 лет, Elderly - добровольцы в возрасте 60 лет.


<table style='NAborder-bottom: 0; font-size: 14px; font-family: "Source Sans Pro", helvetica, sans-serif; width: auto !important; ' class=" lightable-classic">
<caption style="font-size: initial !important;"><b>Table 1. Baseline characteristics of the study participants.</b></caption>
 <thead>
  <tr>
   <th style="text-align:left;font-weight: bold;"> Characteristic </th>
   <th style="text-align:center;font-weight: bold;"> Low Responder<br>N = 14 </th>
   <th style="text-align:center;font-weight: bold;"> High Responder<br>N = 14 </th>
   <th style="text-align:center;font-weight: bold;"> p-value </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">Study arm, n (%)</span> </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> 0.023 </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Young </td>
   <td style="text-align:center;"> 4 (29%) </td>
   <td style="text-align:center;"> 10 (71%) </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Elderly </td>
   <td style="text-align:center;"> 10 (71%) </td>
   <td style="text-align:center;"> 4 (29%) </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">Gender, n (%)</span> </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> 0.430 </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Female </td>
   <td style="text-align:center;"> 8 (57%) </td>
   <td style="text-align:center;"> 10 (71%) </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Male </td>
   <td style="text-align:center;"> 6 (43%) </td>
   <td style="text-align:center;"> 4 (29%) </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">Race, n (%)</span> </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> 0.516 </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> White </td>
   <td style="text-align:center;"> 11 (79%) </td>
   <td style="text-align:center;"> 8 (57%) </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Black or African American </td>
   <td style="text-align:center;"> 2 (14%) </td>
   <td style="text-align:center;"> 5 (36%) </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Unknown </td>
   <td style="text-align:center;"> 1 (7.1%) </td>
   <td style="text-align:center;"> 1 (7.1%) </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">Ethnicity, n (%)</span> </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> &gt;0.999 </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Not Hispanic or Latino </td>
   <td style="text-align:center;"> 12 (86%) </td>
   <td style="text-align:center;"> 13 (93%) </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Hispanic or Latino </td>
   <td style="text-align:center;"> 2 (14%) </td>
   <td style="text-align:center;"> 1 (7.1%) </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">log2(IgG), ELISA</span> </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> 0.012 </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Mean (SD) </td>
   <td style="text-align:center;"> 14.2 (1.0) </td>
   <td style="text-align:center;"> 13.0 (1.3) </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Median (Q1-Q3) </td>
   <td style="text-align:center;"> 14.0 (13.4-15.1) </td>
   <td style="text-align:center;"> 13.0 (12.6-13.7) </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Range </td>
   <td style="text-align:center;"> 12.8-15.9 </td>
   <td style="text-align:center;"> 10.6-15.4 </td>
   <td style="text-align:center;">  </td>
  </tr>
</tbody>
<tfoot><tr><td style="padding: 0; " colspan="100%">
<sup>1</sup> p-value: Study arm, Gender, Race - Pearson's Chi-squared test;<br> 
Ethnicity - Fisher's exact test; IgG - Mann-Whitney test</td></tr></tfoot>
</table>

<br>

## **Описательная статистика: экспрессии генов**

<br>

Исходный датасет включал в себя данные по экспрессии 26925 генов, для нашего исследования удалим те гены, данные по экспрессии которых были пропущены для всех участников исследования, включенных в анализ, - останется матрица экспрессий по 16146 генам. По данным об экспрессии этих генов во _всех_ точках исследования оценим медианное абсолютное отклонение (MAD) и оставим для дальнейшего анализа 5 тыс. генов с наибольшим его значением (гены с наибольшей вариацией экспрессии).

<br>

Ниже представлены гистограммы для экспрессий во всех точках исследования для 10 генов, случайно выбранных из оставшихся 5 тысяч. 

_Если вы запускаете этот отчёт в RStudio и активизируете в опциях .rmd-файла (yaml) опцию `runtime: shiny`, то сможете генерировать случайные выборки по 10 генов по нажатию кнопки New sample и видеть соответствующие графики_

![](ImmunoSpace_SDY984_files/figure-html/fig1_hists-1.png)<!-- -->

<br>

## **Сравнение экспрессии генов между группами испытуемых с разным уровнем ответа: попарные тесты Манна-Уитни**

<br>

В каждой точке исследования по каждому гену сравним экспрессии генов между группами со слабым и сильным ответом на вакцинацию с помощью теста Манна-Уитни, скорректируем полученные p-values по методу Бенджамини-Хохберга (для контроля FDR, в каждой точке отдельно) и посмотрим, по какому количеству генов средние значения экспрессий различаются на уровне значимости 5% без коррекции и с коррекцией.

<table class=" lightable-classic" style='font-size: 14px; font-family: "Source Sans Pro", helvetica, sans-serif; width: auto !important; '>
<caption style="font-size: initial !important;">
<b>Table 2. P-values of the Mann-Whitney tests before and after<br>
                 adjustment under Benjamini &amp; Hochberg method</b>,<br>number of genes (%)</caption>
 <thead>
  <tr>
   <th style="text-align:left;font-weight: bold;">  </th>
   <th style="text-align:center;font-weight: bold;"> Baseline </th>
   <th style="text-align:center;font-weight: bold;"> 1 d. </th>
   <th style="text-align:center;font-weight: bold;"> 3 d. </th>
   <th style="text-align:center;font-weight: bold;"> 7 d. </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">P-value</span> </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> &lt; 0.05 </td>
   <td style="text-align:center;"> 356 (7.1%) </td>
   <td style="text-align:center;"> 423 (8.5%) </td>
   <td style="text-align:center;"> 353 (7.1%) </td>
   <td style="text-align:center;"> 426 (8.5%) </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> $\geq$ 0.05 </td>
   <td style="text-align:center;"> 4644 (93%) </td>
   <td style="text-align:center;"> 4577 (92%) </td>
   <td style="text-align:center;"> 4647 (93%) </td>
   <td style="text-align:center;"> 4574 (91%) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">Adjusted p-value</span> </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> &lt; 0.05 </td>
   <td style="text-align:center;"> 0 (0%) </td>
   <td style="text-align:center;"> 0 (0%) </td>
   <td style="text-align:center;"> 0 (0%) </td>
   <td style="text-align:center;"> 2 (&lt;0.1%) </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> $\geq$ 0.05 </td>
   <td style="text-align:center;"> 5000 (100%) </td>
   <td style="text-align:center;"> 5000 (100%) </td>
   <td style="text-align:center;"> 5000 (100%) </td>
   <td style="text-align:center;"> 4998 (100%) </td>
  </tr>
</tbody>
</table>

<br>

Оказалось всего 2 гена, экспрессии которых статистически значимым образом (на 5%-ном уровне значимости) различаются между сильными и слабыми респондерами, оба наблюдения относятся к точке 7 дней после вакцинации, - это гены IL2RA, RIC3. Статистика по экспрессии этих генов в группах в каждой точке исследования представлена в таблице ниже.

<table style='NAborder-bottom: 0; font-size: 11px; font-family: "Source Sans Pro", helvetica, sans-serif; width: auto !important; ' class=" lightable-classic">
<caption style="font-size: initial !important;"><b>Table 3. Descriptive statistics for genes with expressions,
  significantly different between low and high responders.</b></caption>
 <thead>
<tr>
<th style="empty-cells: hide;border-bottom:hidden;" colspan="1"></th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="3"><div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; "><b>Baseline</b></div></th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="3"><div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; "><b>1 d.</b></div></th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="3"><div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; "><b>3 d.</b></div></th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="3"><div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; "><b>7 d.</b></div></th>
</tr>
  <tr>
   <th style="text-align:left;font-weight: bold;"> Gene/ Statistic </th>
   <th style="text-align:center;font-weight: bold;"> Low Responder<br>N = 14 </th>
   <th style="text-align:center;font-weight: bold;"> High Responder<br>N = 14 </th>
   <th style="text-align:center;font-weight: bold;"> p-value </th>
   <th style="text-align:center;font-weight: bold;"> Low Responder<br>N = 14 </th>
   <th style="text-align:center;font-weight: bold;"> High Responder<br>N = 14 </th>
   <th style="text-align:center;font-weight: bold;"> p-value </th>
   <th style="text-align:center;font-weight: bold;"> Low Responder<br>N = 14 </th>
   <th style="text-align:center;font-weight: bold;"> High Responder<br>N = 14 </th>
   <th style="text-align:center;font-weight: bold;"> p-value </th>
   <th style="text-align:center;font-weight: bold;"> Low Responder<br>N = 14 </th>
   <th style="text-align:center;font-weight: bold;"> High Responder<br>N = 14 </th>
   <th style="text-align:center;font-weight: bold;"> p-value </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">IL2RA</span> </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> 0.914 </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> 0.902 </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> 0.893 </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> 0.017 </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Mean (SD) </td>
   <td style="text-align:center;"> 4.9 (0.4) </td>
   <td style="text-align:center;"> 5.1 (0.4) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> 5.0 (0.4) </td>
   <td style="text-align:center;"> 5.1 (0.4) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> 5.1 (0.4) </td>
   <td style="text-align:center;"> 5.0 (0.5) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> 4.8 (0.3) </td>
   <td style="text-align:center;"> 5.4 (0.4) </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Median (Q1-Q3) </td>
   <td style="text-align:center;"> 5.1 (4.7-5.1) </td>
   <td style="text-align:center;"> 4.9 (4.8-5.4) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> 5.1 (4.7-5.2) </td>
   <td style="text-align:center;"> 5.2 (4.8-5.4) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> 5.1 (4.8-5.4) </td>
   <td style="text-align:center;"> 4.9 (4.6-5.3) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> 4.8 (4.7-5.0) </td>
   <td style="text-align:center;"> 5.3 (5.3-5.5) </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Range </td>
   <td style="text-align:center;"> 4.3-5.3 </td>
   <td style="text-align:center;"> 4.5-5.9 </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> 4.3-5.6 </td>
   <td style="text-align:center;"> 4.7-6.0 </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> 4.3-5.8 </td>
   <td style="text-align:center;"> 4.4-5.8 </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> 4.1-5.2 </td>
   <td style="text-align:center;"> 4.8-6.2 </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">RIC3</span> </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> 0.752 </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> 0.401 </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> 0.739 </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> 0.017 </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Mean (SD) </td>
   <td style="text-align:center;"> 4.2 (0.3) </td>
   <td style="text-align:center;"> 4.4 (0.3) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> 4.1 (0.3) </td>
   <td style="text-align:center;"> 4.4 (0.5) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> 4.2 (0.5) </td>
   <td style="text-align:center;"> 4.5 (0.4) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> 3.8 (0.3) </td>
   <td style="text-align:center;"> 4.5 (0.3) </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Median (Q1-Q3) </td>
   <td style="text-align:center;"> 4.1 (4.0-4.4) </td>
   <td style="text-align:center;"> 4.3 (4.2-4.6) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> 4.1 (3.9-4.2) </td>
   <td style="text-align:center;"> 4.4 (4.2-4.7) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> 4.1 (3.8-4.6) </td>
   <td style="text-align:center;"> 4.5 (4.3-4.6) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> 3.7 (3.6-3.9) </td>
   <td style="text-align:center;"> 4.4 (4.2-4.7) </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Range </td>
   <td style="text-align:center;"> 3.7-4.7 </td>
   <td style="text-align:center;"> 3.8-5.0 </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> 3.5-4.6 </td>
   <td style="text-align:center;"> 3.1-5.1 </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> 3.4-5.0 </td>
   <td style="text-align:center;"> 3.7-5.2 </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> 3.3-4.4 </td>
   <td style="text-align:center;"> 4.1-4.9 </td>
   <td style="text-align:center;">  </td>
  </tr>
</tbody>
<tfoot><tr><td style="padding: 0; " colspan="100%">
<sup>1</sup> p-value: Mann-Whitney test with Benjamini &amp; Hochberg correction</td></tr></tfoot>
</table>

<br>

В качестве иллюстрации:

![](ImmunoSpace_SDY984_files/figure-html/fig2_geneexpr-1.png)<!-- -->

<br>

[**IL2RA**](https://en.wikipedia.org/wiki/IL2RA) - это interleukin 2 receptor subunit alpha ген, который кодирует белок Interleukin-2 receptor alpha chain (CD25), [Interleukin-2 (IL-2)](https://en.wikipedia.org/wiki/Interleukin_2), участвующий в регулировании активности лейкоцитов, отвечающих за иммунитет. Все биологические процессы, в которых участвуют белковые продукты данного гена (по данным из Gene Ontology): inflammatory response to antigenic stimulus; regulation of T cell tolerance induction; apoptotic process; activation-induced cell death of T cells; inflammatory response; immune response; cell surface receptor signaling pathway; Notch signaling pathway; interleukin-2-mediated signaling pathway; positive regulation of activated T cell proliferation; negative regulation of T cell proliferation; positive regulation of T cell differentiation; regulation of T cell homeostatic proliferation; negative regulation of inflammatory response.

[**RIC3**](https://en.wikipedia.org/wiki/RIC3) - это RIC3 acetylcholine receptor chaperone ген, который кодирует chaperon белок RIC-3, отвечающий за резистентность к ингибиторам cholinesterase 3. [Chaperon белки](https://en.wikipedia.org/wiki/Chaperone_(protein)) участвуют в сворачивании или разворачивании крупных белков или макромолекулярных белковых комплексов. Все биологические процессы, в которых участвуют белковые продукты данного гена (по данным из Gene Ontology): protein folding; positive regulation of cytosolic calcium ion concentration; synaptic transmission, cholinergic; protein localization to cell surface; cellular protein-containing complex assembly; positive regulation of protein localization to cell surface.

<br>

## **Сравнение экспрессии генов между группами испытуемых с разным уровнем ответа: линейная смешанная модель**

<br>

### **Описание модели**

<br>

Для каждого гена из ранее отобранных 5 тыс. оценим **линейные смешанные модели (linear mixed effects model) вида**:

1. $expr_{it} = \beta_0 + \beta_{0i} + \beta_1*time_{j} + \beta_2*response_{i} + \epsilon_{ij}$, где:

- $expr_{ij}$ - экспрессия гена для _i_-го испытуемого в точке исследования _j_; оценим отдельно спецификации с:

  - исходными данными по экспрессии, 

  - рангами генов по экспрессии, рассчитанными для каждого испытуемого в каждой точке исследования как квантиль эмпрической функции распределения экспрессий по всем генам для этого испытуемого в этой точке, 

- $time_j$ - точка исследования (в днях от вакцинации), 

- $response_i$ - ответ на вакцинацию _i_-го испытумого (1 - high responder, 0 - low responder), 

- $\beta_0$ - константа регрессионного уравнения (глобальное среднее значение экспрессии гена по всем испытуемым во всех точках),

- $\beta_{0i}$ - случайный эффект (моделируется для каждого испытуемого, имеет нормальное распределение со средним 0 и дисперсией, которая характеризует вариацию средней экспрессии гена между испытуемыми вокруг глоабального среднего значения),

- $\beta_{1,2}$ - коэффициенты регрессионного уравнения (фиксированные эффекты времени от вакцинации и ответа на неё), 

- $\epsilon_{ij}$ - случайный остаток регрессионного уравнения (предполагается, что имеет нормальное распределение со средним значением 0 и дисперсией, которая характеризует вариацию экспрессии для каждого испытуемого (во времени)).

Оценка коэффициента $\beta_2$ данного регрессионного уравнения может интерпретироваться как **средний предельный эффект (average marginal/ treatment effect, AME/ ATE)** ответа на вакцинацию в отношении экспрессии, т.е. в нашем случае она покажет, на сколько отличается среднее значение экспрессии гена среди high респондеров по сравнению с low респондерами по всем точкам исследования.

<br>

2. $expr_{it} = \beta_0 + \beta_{0i} + \beta_1*time_{j} + \beta_2*response_{i} + \beta_3*time_j * response_i + \epsilon_{ij}$, где по сравнению с предыдущей спецификацией добавляется $time_j * response_i$ - пересечение времени анализа и статуса по ответу на вакцинацию.

При оценке этой модели предполагается, что эффект ответа на вацинацию может меняться во времени (и, аналогично, динамика экспрессии может быть разной в зависимости от ответа на вакцинацию). Оценка коэффициента $\beta_2$ в этой модели может быть интерпретирована как ATE в точке baseline (когда $time=0$), т.е. как baseline разница в средней экспрессии гена между теми, кто в дальнейшем был отнесен к high респондерам, и теми, кто нет. Для остальных временных точек ATE может быть получен из соотношения $\beta_2 + \beta_3*time$. Соответственно оценку коэффициента $\beta_3$ мы можем интерпретировать как изменение ATE ответа при увеличении времени от вакцинации на 1 день.

<br>

Оценка линейных смешанных моделей будет выполнена с помощью функции `lmer` из пакета `lme4`, а оценка средних предельных эффектов - с помощью функции `marginaleffects` из одноименного пакета. Заметим, что если в `marginaleffects` после модели с эффектами пересечения задать конкретные временные точки (например, $time=$ 0,1,3,7), то для переменной $response$ мы получим искомые значения ATE в этих точках ($\beta_2 + \beta_3*time$) вместе с их p-values, а если этого не сделать, то на выходе для переменной $response$ мы получим оценку ATE, равную оценке коэффициента $\beta_2$ из модели без эффекта пересечения.

<br>

### **Результаты оценки**

<br>

#### **Наборы значимых генов в разных моделях**

<br>

Всего получилось 414 генов в модели с исходными значениями для экспрессии и 386 генов в модели с рангами для экспрессии, для которых на 5%-ном уровне значимости статистически значим коэффициент при ответе **в регрессиях без эффекта пересечения** - иными словами, значимо различается среднее значение экспрессии между испытуемыми с разным уровнем ответа по всем точкам исследования. 

**В регрессиях с эффектами пересечения** обнаружилось 49 генов в модели с исходными значениями для экспрессии и 42 гена в модели с рангами для экспрессии, для которых на 5%-ном уровне значимости статистически значим и основной эффект для ответа, и эффект пересечения ответа со временем, т.е. испытуемые с разным уровнем ответа значимо отличались не только по среднему значению экспрессии данных генов до вакцинации, но и по его динамике. Ещё для 277 (260) генов в соответствующих моделях был значим только основной эффект, т.е. испытуемые с разным ответом на вакцинацию изначально различались по экспрессии этих генов, но её динамика была схожей в этих группах. Ещё для 318 (337) генов был значим только эффект пересечения, т.е. испытуемые с разным ответом на вакцинацию изначально не различались по экспрессии этих генов, но её динамика была в этих группах разной.

Мы можем **сопоставить p-values при основном эффекте и эффекте пересечения** между ответом и временем анализа для каждого гена - например, с помощью диаграммы рассеяния: по оси X покажем p-value для основного эффекта ответа (оценки коэффициента $\beta_2$), по оси Y - для эффекта пересечения ответа со временем (оценки коэффициента $\beta_3$). Для большей наглядности будем использовать $-log_{10}$-преобразование для обоих p-values. В каждом квадранте, за исключением нижнего левого (сюда попадают гены, у которых оба коэффициента были статистически незначимыми), подпишем по 5 генов с наименьшими значениями из обоих p-values.

![](ImmunoSpace_SDY984_files/figure-html/fig3_p12-1.png)<!-- -->

Также мы можем посмотреть, **насколько пересекаются наборы генов-находок**, выявленных по статистической значимости коэффициентов в моделях с и без эффектов пересечения.

![](ImmunoSpace_SDY984_files/figure-html/fig4_upset_int_noint-1.png)<!-- -->

<br>

Из всех уникальных генов, для которых был статистически значим хотя бы один коэффициент в линейных смешанных моделях, **больше всего генов имели значимый коэффициент при эффекте пересечения ответа со временем** (независимо от способа представления данных по экспрессии). 

Всего по результатам обоих подходов к спецификации экспрессии генов в линейных смешанных моделях нашлось **334 гена, для которых был значим основной эффект в одной из этих спецификаций и незначим эффект пересечения ни в одной из них**, а также **370 генов, для которых, наоборот, был значим эффект пересечения в одной из моделей, но не значим основной эффект ни в одной из них**. Можно считать, что средняя экспрессия генов первой группы значимо различалась между испытуемыми с разным уровнем ответа на вакцинацию до её проведения, но её динамика была для них схожей, а для второй группы генов, наоборот, изменение средней экспрессии после вакцинации значимо различалось в зависимости от ответа, хотя до вакцинации различий в средней экспрессии не было.

С помощью базы Gene Ontology для каждого из этих наборов генов мы определили, в какие **сигнальные пути (биологические процессы)** входят их продукты, взяли пути с наибольшим числом генов, их представляющих, и на графике ниже представили результаты, дополнительно выделив те процессы, которые непосредственно связаны с иммунным ответом (на графике они обозначены как related to immune response).

![](ImmunoSpace_SDY984_files/figure-html/fig5_bps_b2b3-1.png)<!-- -->

![](ImmunoSpace_SDY984_files/figure-html/fig5_bps_b2b3-2.png)<!-- -->

<br>

Для моделей с эффектами пересечения мы также можем оценить **средний эффект воздействия (average treatment effect, ATE) ответа в отдельных точках исследования**. Ниже на графике покажем, насколько пересекаются наборы генов-находок между точками исследования (гены с p-value < 0.05 для соответствующего ATE).

![](ImmunoSpace_SDY984_files/figure-html/fig6_upset_ame-1.png)<!-- -->

![](ImmunoSpace_SDY984_files/figure-html/fig6_upset_ame-2.png)<!-- -->

Для относительно небольшого числа генов средний эффект воздействия ответа в отношении экспрессии является статистически значимым на всём протяжении исследования (от 0 до 7 дней). Также относительно невелико число генов, для которых ATE был статистически незначим до вакцинации, но был значимым с 1 по 7 дни, как и генов, для которых, наоборот значимость ATE до вакцинации сменилась на статистическую незначимость с 1 по 7 дней после неё. Но самое примечательное - 2 моды: **для наибольшего числа генов ATE стал значимым только к 7 дням после вакцинации, для следующего по количеству генов набора ATE был значим в 0-3 дня, но стал незначимым в 7 дней**.

Объедними результаты по моделям с исходными данными по экспрессии и с рангами генов по ней - получится, что для 310 генов ATE в 7 дней был значим в одной из этих моделей, но не значим ни в одной из них в 0-3 дня. Для 216 генов ATE был значим одновременно в 0-1-3 дня в любой из этих моделей и незначим в 7 дней ни в одной из них. Для каждого из этих наборов генов отберём наиболее представленные **биологические процессы**, в которых участвуют их продукты, и покажем их на графике ниже.

![](ImmunoSpace_SDY984_files/figure-html/fig7_bps_ame-1.png)<!-- -->

![](ImmunoSpace_SDY984_files/figure-html/fig7_bps_ame-2.png)<!-- -->

<br>

#### **Статистическая значимость vs. размер эффекта**

<br>

Следующим шагом в анализе результатов оценки линейных смешанных моделей будет **сопоставление статистической значимости среднего эффекта ответа на вакцинацию с размером этого эффекта**. 

Ниже представлен volcano plot для результатов оценки регрессий без эффектов пересечения, где по оси Y показано значение $-log_{10}$-преобразования p-value, а по оси X - оценка соответствующего коэффициента при переменной для ответа. Оранжевым выделим по 5 генов с p-value < 0.05 и самыми низкими и самыми высокими значениями оценки коэффициента.

![](ImmunoSpace_SDY984_files/figure-html/fig8_volcano-1.png)<!-- -->

<br>

Аналогичные графики можно построить и для оценки **среднего эффекта ответа в отдельных точках исследования**, используя результаты для регрессий с эффектами пересечения. На каждом графике оранжевым выделим по 5 генов с p-value < 0.05 и самыми низкими и самыми высокими значениями оценки ATE.

![](ImmunoSpace_SDY984_files/figure-html/fig9_volcano_time-1.png)<!-- -->

<br>

С одной стороны, визуально заметно некоторое "растягивание" точек вправо и влево к 7 дням от вакцинации, а с другой стороны, по отдельным генам, для которых разница в средней экспрессии до вакцинации между сильными и слабыми респондервами была наибольшей по абсолютной величине, заметно её уменьшение.

<br>

С помощью базы Gene Ontology посмотрим, **в какие биологические процессы входят продукты генов, для которых ATE оказался статистически значимым и большим по абсолютной величине**, при этом назовем upregulated гены, для которых ATE > 0.25 для исходных данных по экспрессии генов, downregulated - для которых ATE < -0.25. Результаты для регрессий с рангами по экспрессии здесь учитывать не будем. Покажем на графиках ниже результаты, полученные после отбора по 10 процессов с наибольшим числом генов в каждой точке исследования, выделив отдельно процессы, являющиеся потомками процессов immune system process или defense response.

![](ImmunoSpace_SDY984_files/figure-html/fig10_bps_big-1.png)<!-- -->

<br>

А теперь отдельно только **для процессов, непосредственно связанных с иммунным ответом** (с предварительным отбором до 15 наиболее представленных процессов в каждой точке):

![](ImmunoSpace_SDY984_files/figure-html/fig11_bps_big_imm-1.png)<!-- -->

<br>

## **Взаимосвязь между вероятностью ответа и экспрессией генов: логистическая регрессия**

<br>

### **Описание модели**

<br>

Для каждого гена из ранее отобранных 5 тыс. оценим **логистическую регрессию вида**:

$log(odds(response_{i})) = \gamma_0 + \gamma_1*expr_{0i} + \gamma_2*expr_{1i} + \gamma_3*expr_{3i} + \gamma_4*expr_{7i}$, где:

- $response_i$ - ответ на вакцинацию _i_-го испытумого (1 - high responder, 0 - low responder), 

- $expr_{ji}$ - экспрессия гена для _i_-го испытуемого в точке исследования _j_, $j=0,1,3,7$. Как и при оценке дифференциальной экспрессии, отдельно оценим спецификации с:

  - исходными данными по экспрессии, 

  - рангами генов по экспрессии, рассчитанными для каждого испытуемого в каждой точке исследования как квантиль эмпрической функции распределения экспрессий по всем генам для этого испытуемого в этой точке (в данном случае для удобства интерпретации - см. ниже - мы переведём ранги в шкалу от 0 до 100), 

- $\gamma_0$ - константа регрессионного уравнения,

- $\gamma_1,\gamma_2,\gamma_3,\gamma_4$ - коэффициенты регрессионного уравнения.

Результаты оценки логистических регрессий будем представлять в виде экспонированных значений коэффициентов, которые можно будет интерпретировать как количество раз, в которое изменится шанс сильного ответа на вакцинацию в случае увеличения экспрессии в соответствующей точке на 1, при прочих равных условиях. 

<br>

### **Результаты оценки**

<br>

#### **Наборы значимых генов в разных точках**

<br>

Ниже на графиках покажем, **насколько пересекаются наборы генов-находок**, выявленных по статистической значимости коэффициентов при переменных для экспрессии в разных точках исследования в логистических регрессиях.

![](ImmunoSpace_SDY984_files/figure-html/fig12_upset_logreg-1.png)<!-- -->

![](ImmunoSpace_SDY984_files/figure-html/fig12_upset_logreg-2.png)<!-- -->

Также может быть интересно **сопоставить для каждой точки исследования набор генов со значимым средним эффектом ответа на экспрессию с набором генов со значимым отношением шансов ответа при увеличении экспрессии гена на 1** - результаты на графиках ниже:

![](ImmunoSpace_SDY984_files/figure-html/fig13_compmodels-1.png)<!-- -->

![](ImmunoSpace_SDY984_files/figure-html/fig13_compmodels-2.png)<!-- -->

![](ImmunoSpace_SDY984_files/figure-html/fig13_compmodels-3.png)<!-- -->

![](ImmunoSpace_SDY984_files/figure-html/fig13_compmodels-4.png)<!-- -->

Видим, что наборы, найденные с помощью двух разных подходов, относительно слабо пересекаются друг с другом до вакцинации и в 1-3 дня после неё - и только в 7 дней пересечение становится заметным.

<br>

#### **Статистическая значимость vs. размер эффекта**

<br>

Размер эффекта экспрессии гена в отношении вероятности сильного ответа на вакцинацию можно оценить с помощью **отношения шансов** (экспонированного значения коэффициента логистической регрессии при переменной для экспрессии в соответствующей точке исследования).

Ниже на volcano plot **сопоставим статистическую значимость эффекта с его размером**, для чего по оси Y покажем значение $-log_{10}$-преобразования p-value, а по оси X $log_{10}$ оценки отношения шансов сильного ответа при увеличении соответствующей переменной для экспрессии на 1 (для удобства отображения заменим $log_{10}$ отношения шансов в регрессиях с исходными данными по экспрессии за пределами интервала [-10, 10] на ближайшую из границ этого интервала, в регрессиях для рангов по экспрессии сделаем то же самое по интервалу [-4, 4] - на графиках видно, что все гены, на отображение которых это повлияло, имели статистически незначимые коэффициенты при соответствующих переменных). Оранжевым выделим по 5 генов с p-value < 0.05 и самыми низкими и самыми высокими значениями отношений шансов.

![](ImmunoSpace_SDY984_files/figure-html/fig14_volcano_or-1.png)<!-- -->

<br>

## **Функциональный анализ**

<br>

Выше мы уже использовали базу Gene Ontology (GO) для определения аннотаций биологических процессов, в которых участвуют белковые продукты генов, дифференциально экспрессированных у испытуемых с разным уровнем ответа на вакцинацию. Однако поскольку генов-находок может быть много, искать подобные аннотации по каждому из них не только не слишком удобно, но не совсем правильно, поскольку продукты разных генов могут участвовать в одних и тех же процессах, как и продукты одного и того же гена могут участвовать в разных. Методы, используемые в рамках **singular enrichment анализа (SEA)** (Tipney and Hunter 2010), позволяют определить, какие именно биологические процессы "перепредставлены" (overrepresented) в найденном наборе из $n$ генов по сравнению с полным набором $N$ генов, отобранных нами для анализа (в нашем случае 5 тыс.). 

В своём исследовании мы будем выявлять такие процессы, основываясь на оценке **p-values гипергеометрического распределения**, т.е. вероятностей получить в выборке из $n$ генов такую же или бОльшую представленность данного процесса, как в полном наборе из $N$ генов. Рассчитать p-value можно по формуле (Boyle et al. 2004): $$P=1-\sum_{i=0}^{k-1}\frac{C_M^i C_{N-M}^{n-i}}{C_N^i}$$
где $C_M^i$ - число сочетаний из $M$ по $i$, $N$ - общее количество генов, отобранных нами для анализа (5 тыс.), $M$ - кол-во генов из числа $N$ с определённым биологическим процессом в аннотации, $n$ - кол-во "значимых" генов, обнаруженных на этапе оценки дифференциальной экспрессии, $k$ - кол-во генов из числа $n$ с определённым биологическим процессом в аннотации. 

В R для этого будем использовать функцию `hyperGtest` из пакета `Category` с опцией `conditional = TRUE` (эта опция означает, что при тестировании родительских "узлов" биологических процессов будут исключаться те "потомки", для которых p-value уже оказался ниже порогового уровня значимости; в качестве такого порога будем использовать значение 0.05). Данная функция - это аналог такой же функции из пакета `GOstats`, позволяющий дополнительно вводить ограничения на минимальный и максимальный размеры наборов генов, участвующих в оценке условного гипергеометрического теста (возьмём 10 и 500, соответственно). 

После оценки гипергеометрического теста исключим выявленные в нём биологические процессы, в которых вошло менее 5 генов из числа значимых - далее будем работать с p-values по оставшимся процессам.

Для **контроля FDR** мы попробуем использовать два подхода: коррекцию p-values, полученных по всем выявленным биологическим процессам, с помощью метода Бенджамина-Хохберга и оценку q-values (Storey and Tibshirani 2003).

SEA будем проводить отдельно для каждого набора значимых генов, выявленных нами по результатам оценки различных линейных смешанных моделей и логистических регрессий.

Ниже для каждого набора значимых генов, выявленного нами по результатам оценки линейных смешанных моделей и логистических регрессий, покажем **количество "перепредставленных" этим набором биологических процессов при разных пороговых значениях для скорректированных p-values и для q-values** (пустая ячейка означает отсутствие процессов, удовлетворяющих соответствующему критерию).

<table class=" lightable-classic" style='font-size: 14px; font-family: "Source Sans Pro", helvetica, sans-serif; width: auto !important; '>
<caption style="font-size: initial !important;"><b>Table 5. Number of significant biological processes (GO) found in hypergeometric tests</b></caption>
 <thead>
<tr>
<th style="empty-cells: hide;border-bottom:hidden;" colspan="1"></th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="4"><div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">Adjusted p-value</div></th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="4"><div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">q-value</div></th>
</tr>
  <tr>
   <th style="text-align:left;font-weight: bold;"> Gene set source </th>
   <th style="text-align:center;font-weight: bold;"> &lt;  0.1 </th>
   <th style="text-align:center;font-weight: bold;"> &lt;  0.05 </th>
   <th style="text-align:center;font-weight: bold;"> &lt;  0.01 </th>
   <th style="text-align:center;font-weight: bold;"> &lt;  0.001 </th>
   <th style="text-align:center;font-weight: bold;"> &lt;  0.1 </th>
   <th style="text-align:center;font-weight: bold;"> &lt;  0.05 </th>
   <th style="text-align:center;font-weight: bold;"> &lt;  0.01 </th>
   <th style="text-align:center;font-weight: bold;"> &lt;  0.001 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Expr., sig. b2 (Baseline) </td>
   <td style="text-align:center;"> 8 </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> 20 </td>
   <td style="text-align:center;"> 12 </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Expr., sig. b3 (Change) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Expr., sig. b2 (Overall) </td>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Expr., sig. ATE (Baseline) </td>
   <td style="text-align:center;"> 3 </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> 10 </td>
   <td style="text-align:center;"> 4 </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Expr., sig. ATE (1 d.) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Expr., sig. ATE (3 d.) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Expr., sig. ATE (7 d.) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Ranks, sig. b2 (Baseline) </td>
   <td style="text-align:center;"> 6 </td>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> 55 </td>
   <td style="text-align:center;"> 7 </td>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Ranks, sig. b3 (Change) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Ranks, sig. b2 (Overall) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Ranks, sig. ATE (Baseline) </td>
   <td style="text-align:center;"> 3 </td>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> 39 </td>
   <td style="text-align:center;"> 8 </td>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Ranks, sig. ATE (1 d.) </td>
   <td style="text-align:center;"> 10 </td>
   <td style="text-align:center;"> 8 </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> 38 </td>
   <td style="text-align:center;"> 10 </td>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Ranks, sig. ATE (3 d.) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Ranks, sig. ATE (7 d.) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Expr., sig.OR (Baseline) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> 18 </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Expr., sig.OR (1 d.) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> 136 </td>
   <td style="text-align:center;"> 21 </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Expr., sig.OR (3 d.) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Expr., sig.OR (7 d.) </td>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Ranks, sig.OR (Baseline) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Ranks, sig.OR (1 d.) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Ranks, sig.OR (3 d.) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Ranks, sig.OR (7 d.) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
  </tr>
</tbody>
</table>

Видим, что, в основном, перепредставленные процессы обнаруживаются до вакцинации или в 1 день после неё.

В таблице ниже перечислим все **процессы с q-value < 0.05**, сгруппировав их по точкам, в которых они были перепредставлены по результатам хотя бы одной модели. Среди них нет ни одного, который бы относился к процессам, обозначенным нами как непосредственно связанные с иммунным ответом.

<table class=" lightable-classic" style='font-size: 14px; font-family: "Source Sans Pro", helvetica, sans-serif; width: auto !important; '>
<caption style="font-size: initial !important;"><b>Table 6. Overrepresented BPs, by time</b></caption>
 <thead>
  <tr>
   <th style="text-align:left;font-weight: bold;"> Time </th>
   <th style="text-align:center;font-weight: bold;"> GOBPID </th>
   <th style="text-align:left;font-weight: bold;"> BP (GO term) </th>
   <th style="text-align:left;font-weight: bold;"> Definition (GO) </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 1 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0003416 </td>
   <td style="text-align:left;vertical-align:top;"> endochondral bone growth </td>
   <td style="text-align:left;vertical-align:top;"> The increase in size or mass of an endochondral bone that contributes to the shaping of the bone. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 1 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0006636 </td>
   <td style="text-align:left;vertical-align:top;"> unsaturated fatty acid biosynthetic process </td>
   <td style="text-align:left;vertical-align:top;"> The chemical reactions and pathways resulting in the formation of an unsaturated fatty acid, any fatty acid containing one or more double bonds between carbon atoms. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 1 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0007166 </td>
   <td style="text-align:left;vertical-align:top;"> cell surface receptor signaling pathway </td>
   <td style="text-align:left;vertical-align:top;"> A series of molecular signals initiated by activation of a receptor on the surface of a cell. The pathway begins with binding of an extracellular ligand to a cell surface receptor, or for receptors that signal in the absence of a ligand, by ligand-withdrawal or the activity of a constitutively active receptor. The pathway ends with regulation of a downstream cellular process, e.g. transcription. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 1 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0007569 </td>
   <td style="text-align:left;vertical-align:top;"> cell aging </td>
   <td style="text-align:left;vertical-align:top;"> An aging process that has as participant a cell after a cell has stopped dividing. Cell aging may occur when a cell has temporarily stopped dividing through cell cycle arrest (GO:0007050) or when a cell has permanently stopped dividing, in which case it is undergoing cellular senescence (GO:0090398). May precede cell death (GO:0008219) and succeed cell maturation (GO:0048469). </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 1 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0009893 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of metabolic process </td>
   <td style="text-align:left;vertical-align:top;"> Any process that activates or increases the frequency, rate or extent of the chemical reactions and pathways within a cell or an organism. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 1 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0009967 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of signal transduction </td>
   <td style="text-align:left;vertical-align:top;"> Any process that activates or increases the frequency, rate or extent of signal transduction. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 1 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0010171 </td>
   <td style="text-align:left;vertical-align:top;"> body morphogenesis </td>
   <td style="text-align:left;vertical-align:top;"> The process in which the anatomical structures of the soma are generated and organized. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 1 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0010959 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of metal ion transport </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate, or extent of metal ion transport. Metal ion transport is the directed movement of metal ions, any metal ion with an electric charge, into, out of or within a cell, or between cells, by means of some agent such as a transporter or pore. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 1 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0019369 </td>
   <td style="text-align:left;vertical-align:top;"> arachidonic acid metabolic process </td>
   <td style="text-align:left;vertical-align:top;"> The chemical reactions and pathways involving arachidonic acid, a straight chain fatty acid with 20 carbon atoms and four double bonds per molecule. Arachidonic acid is the all-Z-(5,8,11,14)-isomer. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 1 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0030003 </td>
   <td style="text-align:left;vertical-align:top;"> cellular cation homeostasis </td>
   <td style="text-align:left;vertical-align:top;"> Any process involved in the maintenance of an internal steady state of cations at the level of a cell. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 1 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0030198 </td>
   <td style="text-align:left;vertical-align:top;"> extracellular matrix organization </td>
   <td style="text-align:left;vertical-align:top;"> A process that is carried out at the cellular level which results in the assembly, arrangement of constituent parts, or disassembly of an extracellular matrix. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 1 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0030336 </td>
   <td style="text-align:left;vertical-align:top;"> negative regulation of cell migration </td>
   <td style="text-align:left;vertical-align:top;"> Any process that stops, prevents, or reduces the frequency, rate or extent of cell migration. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 1 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0040013 </td>
   <td style="text-align:left;vertical-align:top;"> negative regulation of locomotion </td>
   <td style="text-align:left;vertical-align:top;"> Any process that stops, prevents, or reduces the frequency, rate or extent of locomotion of a cell or organism. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 1 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0042325 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of phosphorylation </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate or extent of addition of phosphate groups into a molecule. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 1 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0045937 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of phosphate metabolic process </td>
   <td style="text-align:left;vertical-align:top;"> Any process that activates or increases the frequency, rate or extent of the chemical reactions and pathways involving phosphates. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 1 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0046456 </td>
   <td style="text-align:left;vertical-align:top;"> icosanoid biosynthetic process </td>
   <td style="text-align:left;vertical-align:top;"> The chemical reactions and pathways resulting in the formation of icosanoids, any of a group of C20 polyunsaturated fatty acids. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 1 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0048518 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of biological process </td>
   <td style="text-align:left;vertical-align:top;"> Any process that activates or increases the frequency, rate or extent of a biological process. Biological processes are regulated by many means; examples include the control of gene expression, protein modification or interaction with a protein or substrate molecule. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 1 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0048522 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of cellular process </td>
   <td style="text-align:left;vertical-align:top;"> Any process that activates or increases the frequency, rate or extent of a cellular process, any of those that are carried out at the cellular level, but are not necessarily restricted to a single cell. For example, cell communication occurs among more than one cell, but occurs at the cellular level. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 1 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0048844 </td>
   <td style="text-align:left;vertical-align:top;"> artery morphogenesis </td>
   <td style="text-align:left;vertical-align:top;"> The process in which the anatomical structures of arterial blood vessels are generated and organized. Arteries are blood vessels that transport blood from the heart to the body and its organs. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 1 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0050801 </td>
   <td style="text-align:left;vertical-align:top;"> ion homeostasis </td>
   <td style="text-align:left;vertical-align:top;"> Any process involved in the maintenance of an internal steady state of ions within an organism or cell. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 1 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0051174 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of phosphorus metabolic process </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate or extent of the chemical reactions and pathways involving phosphorus or compounds containing phosphorus. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 1 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0051271 </td>
   <td style="text-align:left;vertical-align:top;"> negative regulation of cellular component movement </td>
   <td style="text-align:left;vertical-align:top;"> Any process that stops, prevents, or reduces the frequency, rate or extent of the movement of a cellular component. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 1 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0051480 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of cytosolic calcium ion concentration </td>
   <td style="text-align:left;vertical-align:top;"> Any process involved in the maintenance of an internal steady state of calcium ions within the cytosol of a cell or between the cytosol and its surroundings. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 1 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0051928 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of calcium ion transport </td>
   <td style="text-align:left;vertical-align:top;"> Any process that activates or increases the frequency, rate or extent of the directed movement of calcium ions into, out of or within a cell, or between cells, by means of some agent such as a transporter or pore. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 1 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0055010 </td>
   <td style="text-align:left;vertical-align:top;"> ventricular cardiac muscle tissue morphogenesis </td>
   <td style="text-align:left;vertical-align:top;"> The process in which the anatomical structures of cardiac ventricle muscle is generated and organized. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 1 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0055065 </td>
   <td style="text-align:left;vertical-align:top;"> metal ion homeostasis </td>
   <td style="text-align:left;vertical-align:top;"> Any process involved in the maintenance of an internal steady state of metal ions within an organism or cell. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 1 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0060415 </td>
   <td style="text-align:left;vertical-align:top;"> muscle tissue morphogenesis </td>
   <td style="text-align:left;vertical-align:top;"> The process in which the anatomical structures of muscle tissue are generated and organized. Muscle tissue consists of a set of cells that are part of an organ and carry out a contractive function. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 1 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0072507 </td>
   <td style="text-align:left;vertical-align:top;"> divalent inorganic cation homeostasis </td>
   <td style="text-align:left;vertical-align:top;"> Any process involved in the maintenance of an internal steady state of divalent cations within an organism or cell. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 1 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0098868 </td>
   <td style="text-align:left;vertical-align:top;"> bone growth </td>
   <td style="text-align:left;vertical-align:top;"> The increase in size or mass of a bone that contributes to the shaping of that bone. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 1 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:1902531 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of intracellular signal transduction </td>
   <td style="text-align:left;vertical-align:top;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 1 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:2000772 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of cellular senescence </td>
   <td style="text-align:left;vertical-align:top;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0002067 </td>
   <td style="text-align:left;vertical-align:top;"> glandular epithelial cell differentiation </td>
   <td style="text-align:left;vertical-align:top;"> The process in which a relatively unspecialized cell acquires specialized features of a glandular epithelial cell. A glandular epithelial cell is a columnar/cuboidal epithelial cell found in a two dimensional sheet with a free surface exposed to the lumen of a gland. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0003416 </td>
   <td style="text-align:left;vertical-align:top;"> endochondral bone growth </td>
   <td style="text-align:left;vertical-align:top;"> The increase in size or mass of an endochondral bone that contributes to the shaping of the bone. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0007187 </td>
   <td style="text-align:left;vertical-align:top;"> G protein-coupled receptor signaling pathway, coupled to cyclic nucleotide second messenger </td>
   <td style="text-align:left;vertical-align:top;"> The series of molecular signals generated as a consequence of a G protein-coupled receptor binding to its physiological ligand, where the pathway proceeds with activation or inhibition of a nucleotide cyclase activity and a subsequent change in the concentration of a cyclic nucleotide. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0007569 </td>
   <td style="text-align:left;vertical-align:top;"> cell aging </td>
   <td style="text-align:left;vertical-align:top;"> An aging process that has as participant a cell after a cell has stopped dividing. Cell aging may occur when a cell has temporarily stopped dividing through cell cycle arrest (GO:0007050) or when a cell has permanently stopped dividing, in which case it is undergoing cellular senescence (GO:0090398). May precede cell death (GO:0008219) and succeed cell maturation (GO:0048469). </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0007605 </td>
   <td style="text-align:left;vertical-align:top;"> sensory perception of sound </td>
   <td style="text-align:left;vertical-align:top;"> The series of events required for an organism to receive an auditory stimulus, convert it to a molecular signal, and recognize and characterize the signal. Sonic stimuli are detected in the form of vibrations and are processed to form a sound. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0015872 </td>
   <td style="text-align:left;vertical-align:top;"> dopamine transport </td>
   <td style="text-align:left;vertical-align:top;"> The directed movement of dopamine into, out of or within a cell, or between cells, by means of some agent such as a transporter or pore. Dopamine is a catecholamine neurotransmitter and a metabolic precursor of noradrenaline and adrenaline. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0044262 </td>
   <td style="text-align:left;vertical-align:top;"> cellular carbohydrate metabolic process </td>
   <td style="text-align:left;vertical-align:top;"> The chemical reactions and pathways involving carbohydrates, any of a group of organic compounds based of the general formula Cx(H2O)y, as carried out by individual cells. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0046661 </td>
   <td style="text-align:left;vertical-align:top;"> male sex differentiation </td>
   <td style="text-align:left;vertical-align:top;"> The establishment of the sex of a male organism by physical differentiation. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0048568 </td>
   <td style="text-align:left;vertical-align:top;"> embryonic organ development </td>
   <td style="text-align:left;vertical-align:top;"> Development, taking place during the embryonic phase, of a tissue or tissues that work together to perform a specific function or functions. Development pertains to the process whose specific outcome is the progression of a structure over time, from its formation to the mature structure. Organs are commonly observed as visibly distinct structures, but may also exist as loosely associated clusters of cells that work together to perform a specific function or functions. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0051057 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of small GTPase mediated signal transduction </td>
   <td style="text-align:left;vertical-align:top;"> Any process that activates or increases the frequency, rate or extent of small GTPase mediated signal transduction. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0051148 </td>
   <td style="text-align:left;vertical-align:top;"> negative regulation of muscle cell differentiation </td>
   <td style="text-align:left;vertical-align:top;"> Any process that stops, prevents, or reduces the frequency, rate or extent of muscle cell differentiation. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0051154 </td>
   <td style="text-align:left;vertical-align:top;"> negative regulation of striated muscle cell differentiation </td>
   <td style="text-align:left;vertical-align:top;"> Any process that stops, prevents, or reduces the frequency, rate or extent of striated muscle cell differentiation. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0051482 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of cytosolic calcium ion concentration involved in phospholipase C-activating G protein-coupled signaling pathway </td>
   <td style="text-align:left;vertical-align:top;"> Any process that increases the concentration of calcium ions in the cytosol that occurs as part of a PLC-activating G protein-coupled receptor signaling pathway. G-protein-activated PLC hydrolyses phosphatidylinositol-bisphosphate (PIP2) to release diacylglycerol (DAG) and inositol trisphosphate (IP3). IP3 then binds to calcium release channels in the endoplasmic reticulum (ER) to trigger calcium ion release into the cytosol. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0051928 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of calcium ion transport </td>
   <td style="text-align:left;vertical-align:top;"> Any process that activates or increases the frequency, rate or extent of the directed movement of calcium ions into, out of or within a cell, or between cells, by means of some agent such as a transporter or pore. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0060348 </td>
   <td style="text-align:left;vertical-align:top;"> bone development </td>
   <td style="text-align:left;vertical-align:top;"> The process whose specific outcome is the progression of bone over time, from its formation to the mature structure. Bone is the hard skeletal connective tissue consisting of both mineral and cellular components. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0098868 </td>
   <td style="text-align:left;vertical-align:top;"> bone growth </td>
   <td style="text-align:left;vertical-align:top;"> The increase in size or mass of a bone that contributes to the shaping of that bone. </td>
  </tr>
</tbody>
</table>

<br>

## **Список литературы**

<br>


Boyle, E. I., S. Weng, J. Gollub, H. Jin, D. Botstein, J. M. Cherry, and G. Sherlock. 2004. “GO::TermFinder–Open Source Software for Accessing Gene Ontology Information and Finding Significantly Enriched Gene Ontology Terms Associated with a List of Genes.” Bioinformatics 20 (18): 3710–15. https://doi.org/10.1093/bioinformatics/bth456.

Desmet, Christophe J., and Ken J. Ishii. 2012. “Nucleic Acid Sensing at the Interface Between Innate and Adaptive Immunity in Vaccination.” Nature Reviews Immunology 12 (7): 479–91. https://doi.org/10.1038/nri3247.

Fourati, Slim, Lewis E. Tomalin, Matthew P. Mulè, Daniel G. Chawla, Bram Gerritsen, Dmitry Rychkov, Evan Henrich, et al. 2022. “Pan-Vaccine Analysis Reveals Innate Immune Endotypes Predictive of Antibody Responses to Vaccination.” Nature Immunology 23 (12): 1777–87. https://doi.org/10.1038/s41590-022-01329-5.

Li, Shuzhao, Nicole L. Sullivan, Nadine Rouphael, Tianwei Yu, Sophia Banton, Mohan S. Maddur, Megan McCausland, et al. 2017. “Metabolic Phenotypes of Response to Vaccination in Humans.” Cell 169 (5): 862–877.e17. https://doi.org/10.1016/j.cell.2017.04.026.

Pollard, Andrew J., and Else M. Bijker. 2021. “A Guide to Vaccinology: From Basic Principles to New Developments.” Nature Reviews Immunology 21 (2): 83–100. https://doi.org/10.1038/s41577-020-00479-7.

Storey, John D., and Robert Tibshirani. 2003. “Statistical Significance for Genomewide Studies.” Proceedings of the National Academy of Sciences 100 (16): 9440–45. https://doi.org/10.1073/pnas.1530509100.

Tipney, Hannah, and Lawrence Hunter. 2010. “An Introduction to Effective Use of Enrichment Analysis Software.” Human Genomics 4 (3): 202. https://doi.org/10.1186/1479-7364-4-3-202.
