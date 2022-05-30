change1 = 1.1
change2 = 0.1
obs1 = 100
obs2 = 0.1
pred1 = change1*obs1
pred2 = change2*obs2

N = 2

obs = c(obs1, obs2)
preds = c(pred1,pred2)

obs_sq = sqrt(sum(obs^2)/N)
pred_sq = sqrt(sum(preds^2)/N)
res = sqrt(sum((obs-preds)^2)/N)

RSS = sum((obs-preds)^2)/N

pbpkI = ((res/obs_sq)+ (res/pred_sq))/2

new_index = (  sum(((preds-obs)/preds)^2)/N + sum(((preds-obs)/obs)^2)/N ) /2

print(RSS)
print(pbpkI)
print(new_index)