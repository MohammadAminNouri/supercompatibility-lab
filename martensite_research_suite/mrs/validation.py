import numpy as np
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, homogeneity_score, completeness_score, v_measure_score

def partition_metrics(truth,pred):
    truth=np.asarray(truth); pred=np.asarray(pred)
    return {
      "ARI": adjusted_rand_score(truth,pred),
      "NMI": normalized_mutual_info_score(truth,pred),
      "homogeneity": homogeneity_score(truth,pred),
      "completeness": completeness_score(truth,pred),
      "V_measure": v_measure_score(truth,pred),
      "true_parents": len(set(truth.tolist())),
      "reconstructed_parents": len(set(pred.tolist())),
    }

def boundary_metrics(edges,truth,pred):
    yt=[]; yp=[]
    for i,j in edges:
        yt.append(truth[i]!=truth[j]); yp.append(pred[i]!=pred[j])
    yt=np.asarray(yt,bool); yp=np.asarray(yp,bool)
    tp=np.sum(yt&yp); fp=np.sum(~yt&yp); fn=np.sum(yt&~yp)
    p=tp/(tp+fp) if tp+fp else 1.0
    r=tp/(tp+fn) if tp+fn else 1.0
    f=2*p*r/(p+r) if p+r else 0.0
    jac=tp/(tp+fp+fn) if tp+fp+fn else 1.0
    return {"boundary_precision":p,"boundary_recall":r,"boundary_F1":f,"boundary_Jaccard":jac}
