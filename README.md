# env settings

python3에서 동작을 확인하였습니다
자세한 스펙은 `cat Pipfile`을 참고하시기 바랍니다
terminal에서 python3 커맨드를 사용할 수 있다면 ($PROJECT_ROOT)/main.py DNAseq.fasta만으로 실행됩니다


```
Hairpin Finder by soo. Please Enjoy!

positional arguments:
  input_file     input file to read a DNA sequence

optional arguments:
  -h, --help     show this help message and exit
  -v, --verbose  increase output verbosity
```

### FYI

이 프로젝트는 pipenv를 사용합니다 즉 Pipfile을 활용하고 있습니다
다행히도 서드파티 의존성은 없으므로 신경쓰지 않아도 되지만, 서드파티를 추가할 때는 Pipfile을 활용해 주시기 바랍니다

# How to run
`python3 main.py DNAseq.fasta`
