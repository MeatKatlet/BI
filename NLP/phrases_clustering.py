# -*- coding: utf-8 -*-
from gensim.models.doc2vec import  Doc2Vec, TaggedDocument
import pandas as pd
from sklearn.cluster import KMeans
import random
import statsmodels.api as sm
from collections import OrderedDict
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.metrics.pairwise import cosine_similarity
from scipy.stats import chi2_contingency
from nltk.stem.snowball import SnowballStemmer
import re
import nltk
import numpy as np
from sklearn import metrics
from sklearn.metrics import confusion_matrix
from NLP.supervized_clustering_metrics import precison_recall_f1

#todo прочитать json файл

df = pd.read_json("test_task.json", orient='columns',encoding='utf-8')
#todo посчитать в нем количество групп
clusters = df.shape[0]
max_group_size = df.shape[1]
#todo пометить все предложения в кластерах!, создать колонку во фрейме где будут id предложенгий
sentenceLabeled = []
#df[:1]

before = {}
before_sent = {}

documents = []

sentenceID = 0
b = np.array([])
for group in range(0,clusters,1):
    for col in range(0, max_group_size, 1):
        if (df[col][group] != None):
            sentenceL = TaggedDocument(words=df[col][group].lower().split(), tags=['SENT_%s' % sentenceID])
            # sentenceL = TaggedDocument(words=sentence.split(), tags = ['SENT_%s' %sentenceID])
            sentenceLabeled.append(sentenceL)
            before['SENT_%s' % sentenceID] = group
            before_sent['SENT_%s' % sentenceID] = df[col][group]
            documents.append(df[col][group])
            b = np.append(b, int(group))
            sentenceID += 1


def prepare_data(docs):
    documents2 = []
    stemmer = SnowballStemmer("english")
    stemmer2 = SnowballStemmer("russian")
    for doc in docs:
        doc2 = doc.lower()
        tokens = [word for word in nltk.word_tokenize(doc2)]
        filtered_tokens = []
        # filter out any tokens not containing letters (e.g., numeric tokens, raw punctuation)
        for token in tokens:
            if re.search('[a-zа-я]', token):
                filtered_tokens.append(token)
        stems = [stemmer.stem(t) for t in filtered_tokens]
        stems2 = [stemmer2.stem(t) for t in stems]

        res = ' '.join(stems2)
        documents2.append(res)
    return documents2

prep_docs = prepare_data(documents)

"""
max_df=0.5, max_features=200000,
                                 min_df=0.2, stop_words='english',
                                 use_idf=True, ngram_range=(1,3)
"""
tfidf_vectorizer = TfidfVectorizer()
tfidf_matrix = tfidf_vectorizer.fit_transform(prep_docs)
print(tfidf_matrix.shape)

#cosine_similarity(tfidf_matrix[0:1], tfidf_matrix)

#идем по каждому кластеру и по каждому предложению в кластере, помечаем его

#todo перемешать в нем все предложения
"""
new_list = list(sentenceLabeled)
random.shuffle(new_list)
#todo обучить модель
alpha_val = 0.025        # Initial learning rate
min_alpha_val = 1e-4     # Minimum for linear learning rate decay
passes = 15              # Number of passes of one document during training

alpha_delta = (alpha_val - min_alpha_val) / (passes - 1)

model = Doc2Vec( vector_size = 20 # Model initialization
    , window = 20
    , min_count = 5
    , workers = 1)

model.build_vocab(new_list) # Building vocabulary

for epoch in range(passes):

    # Shuffling gets better results

    random.shuffle(new_list)

    # Train

    #model.alpha, model.min_alpha = alpha_val, alpha_val

    model.train(new_list, total_examples=model.corpus_count, epochs=15)

    # Logs

    print('Completed pass %i at alpha %f' % (epoch + 1, alpha_val))

    # Next run alpha

    #alpha_val -= alpha_delta
"""
"""
model = Doc2Vec(alpha=0.025, min_alpha=0.025,dm=1, vector_size=100)  # use fixed learning rate
model.build_vocab(new_list)
for epoch in range(10):
    model.train(new_list, epochs=10, total_examples=model.corpus_count)
    model.alpha -= 0.002  # decrease the learning rate
    model.min_alpha = model.alpha  # fix the learning rate, no decay
#most_similar для предложения из кластера для самоконтроля
"""

#todo kmeans

km = KMeans(n_clusters=clusters)
km.fit(tfidf_matrix)#????????????
clusters = km.labels_.tolist()
#запустить, посмотреть как это выглядит

#indices = [0, 1, 4]
#np.take(a, indices)
#идем по списку из
"""
indeces_of_groups = {}

for i in range(0,len(clusters),1):
    val = clusters[i]
    if val not in indeces_of_groups:
        indeces_of_groups[val] = np.array([])
        indeces_of_groups[val] = np.append(indeces_of_groups[val], int(i))
    else:
        indeces_of_groups[val] = np.append(indeces_of_groups[val], int(i))


predicted = np.array([])
prev = -1
sentenceID = 0


for i in range(0,len(before),1):
    if before['SENT_'+str(i)] != prev:
        l = len(indeces_of_groups[clusters[i]])
        # indices = [0, 1, 4]

        a = np.full((1, l), before['SENT_'+str(i)])
        #predicted.append(a.fill(before['SENT_'+i]))
        predicted = np.append(predicted, a)
        prev = before['SENT_'+str(i)]
        sentenceID += 1
"""

"""
already = {}
u = np.unique(b)
j = 0
for i in range(0,len(clusters),1):
    val = clusters[i]
    if val not in already:
        already[val] = u[j]
        j += 1


for i in range(0,len(clusters),1):
    clusters[i] = already[clusters[i]]


print(confusion_matrix(b, clusters))
"""
print("Homogeneity: %0.3f" % metrics.homogeneity_score(b, clusters))
print("Completeness: %0.3f" % metrics.completeness_score(b, clusters))
print("V-measure: %0.3f" % metrics.v_measure_score(b, clusters))
print("Adjusted Rand-Index: %.3f" % metrics.adjusted_rand_score(b, clusters))

indeces_of_groups = {}
dict_list = []
for i in range(0,len(clusters),1):
    val = clusters[i]
    if val not in indeces_of_groups:
        indeces_of_groups[val] = []
        indeces_of_groups[val].append(str(b[i]))# = #np.append(indeces_of_groups[val], str(b[i]))
    else:
        indeces_of_groups[val].append(str(b[i]))#indeces_of_groups[val] = np.append(indeces_of_groups[val], str(b[i]))

[dict_list.append(v) for k,v in indeces_of_groups.items()]

precison_recall_f1(dict_list) #cluster_list = [['a','a','a','a','a','b'],['a','b','b','b','b','c'],['a','a','c','c','c']]


#todo сравнить все полученные кластеры с изначальнымми
## Print Sentence Clusters ##
"""
cluster_info = {'sentence': sentenceLabeled, 'cluster' : clusters}
sentenceDF = pd.DataFrame(cluster_info, index=[clusters], columns = ['sentence','cluster'])

for num in range(num_clusters):
     print()
     print("Sentence cluster %d: " %int(num+1), end='')
     print()
     for sentence in sentenceDF.ix[num]['sentence'].values.tolist():
        print(' %s ' %sentence, end='')
        print()
    print()
"""
"""
after = {}

for sentence in range(0,sentenceID,1):
    after[new_list[sentence].tags[0]] = clusters[sentence]

after = OrderedDict(sorted(after.items()))
before = OrderedDict(sorted(before.items()))

table = pd.crosstab(before, after, margins=True)

chi2, p, dof, ex = chi2_contingency(table)
a = 1
"""


#найти кластеры - берем кластер - ищем предложения из него в первоначальных кластерах
#хи квадрат для двух кластерных решений

#4.944975233793541e-36

#todo попробовать убрать цифры из предложений, русские стоп слова - предлоги! как убрать?
#







