set +e
for i in {1..5}
do
     python ReadSampler.py -t 22q11del -p m 22q11del$i  
     python ReadSampler.py -t 22q11del -p p 22q11del$i  
     #python ReadSampler.py -t 22q11dup -p m 22q11del$i
     #python ReadSampler.py -t 22q11dup -p p 22q11del$i
     #python ReadSampler.py -t 22q13del -p m 22q11del$i
     #python ReadSampler.py -t 22q13del -p p 22q11del$i
     #python ReadSampler.py -t complete -p m 22q11del$i
     #python ReadSampler.py -t complete -p p 22q11del$i
     #python ReadSampler.py -t longd -p m 22q11del$i
     #python ReadSampler.py -t longd -p p 22q11del$i
done
