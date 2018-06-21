from collections import Counter

def precison_recall_f1(cluster_list):
    #cluster_list = [['a','a','a','a','a','b'],['a','b','b','b','b','c'],['a','a','c','c','c']]
    num_doc = 0
    positives = 0
    negatives = 0
    TP = 0
    FP = 0
    FN = 0
    TN = 0
    c_list = []
    for c in range(0, len(cluster_list)):
        # calculating num_doc count...
        num_doc += len(cluster_list[c])

        # calculating positives...
        positives += (len(cluster_list[c]) * (len(cluster_list[c]) - 1)) / 2

        # calculating TP...
        c = Counter(cluster_list[c])
        c_list.append(c)
        tp_temp = 0
        for k, v in dict(c).items():
            if v > 1:
                tp_temp += (v * (v - 1)) / 2
        TP += tp_temp

    FP = positives - TP
    negatives = ((num_doc * (num_doc - 1)) / 2) - positives
    # Add all the cluster together
    sum = Counter()
    for c in c_list:
        sum += c
    # calculating FN...
    for ct in c_list:
        fn_temp = 0
        for k, v in dict(ct).items():
            fn_temp += v * (sum[k] - v)
        sum -= ct
        FN += fn_temp
    TN = negatives - FN
    print("num_doc is %d " % num_doc)
    print("positives is %d " % positives)
    print("TP is %d " % TP)
    print("FP is %d " % FP)
    print("FN is %d " % FN)
    print("TN is %d " % TN)

    Precision = TP / (TP + FP)
    print("Precision is %.2f " % Precision)

    Recall = TP / (TP + FN)
    print("Recall is %.2f " % Recall)

    F1 = (2 * Recall * Precision) / (Recall + Precision)


    print("F1 is %.2f " % F1)