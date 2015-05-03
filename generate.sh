set +e
for i in {1..10}
do
    nohup python ReadSampler.py -t 22q11del -p m 22q11del$i.txt & > 22q11delm$1.log
    nohup python ReadSampler.py -t 22q11del -p p 22q11del$i.txt & > 22q11delp$1.log
    nohup python ReadSampler.py -t 22q11dup -p m 22q11del$i.txt & > 22q11dupm$1.log
    nohup python ReadSampler.py -t 22q11dup -p p 22q11del$i.txt & > 22q11dupp$1.log
    nohup python ReadSampler.py -t 22q13del -p m 22q11del$i.txt & > 22q13delm$1.log
    nohup python ReadSampler.py -t 22q13del -p p 22q11del$i.txt & > 22q13delp$1.log
    nohup python ReadSampler.py -t complete -p m 22q11del$i.txt & > completem$1.log
    nohup python ReadSampler.py -t complete -p p 22q11del$i.txt & > completep$1.log
    nohup python ReadSampler.py -t longd -p m 22q11del$i.txt & > longdm$1.log
    nohup python ReadSampler.py -t longd -p p 22q11del$i.txt & > longdp$1.log
done
