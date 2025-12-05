# Multi-Cas Analyzer

CRISPR-Cas9 editing 결과를 분석하기 위한 Python 툴입니다. FASTQ 파일로부터 insertion, deletion, substitution, complex edit 패턴을 자동으로 분석하고 Excel 리포트를 생성합니다.

## 주요 기능

- FASTQ 파일에서 CRISPR 편집 패턴 자동 분석
- Wild Type(WT), Insertion, Deletion, Substitution, Complex edit 분류
- 편집 타입별로 정리된 Excel 리포트 생성
- 시퀀스 정렬(alignment) 결과 시각화
- 다중 타겟 분석 지원

## 요구사항

```bash
pip install -r requirements.txt
```

필요한 패키지:
- biopython
- pandas
- openpyxl
- regex

## 사용법

### 기본 실행

```bash
python multi_cas_analyzer.py <fastq_folder> <parameter_file> <edit_range>
```

### 파라미터

- `fastq_folder`: FASTQ 파일이 있는 폴더 경로
- `parameter_file`: 타겟 정보가 담긴 파라미터 파일
- `edit_range`: cleavage site 주변 분석 범위 (bp)

### 예시

```bash
python multi_cas_analyzer.py ./fastq NovaIscB_VEGFA_wRNA1_param.txt 10
```

## 파라미터 파일 형식

파라미터 파일은 3줄 단위로 타겟 정보를 작성합니다:

```
Target_Name
REFERENCE_SEQUENCE
SPACER_SEQUENCE
```

예시:
```
NovaIscB_VEGFA_wRNA1
ATCGATCGATCG...
ATCGATCG
```

## 출력 결과

프로그램은 `output` 폴더에 다음 파일들을 생성합니다:

1. **텍스트 리포트** (`*.fastq.txt`)
   - 각 FASTQ 파일별 전체 통계
   - WT, Insertion, Deletion, Substitution, Complex 카운트 및 비율

2. **Excel 리포트** (`*_patterns.xlsx`)
   - Summary: 전체 편집 타입 요약
   - 각 편집 타입별 시트
   - Edit range alignment 및 패턴별 카운트

3. **요약 파일** (`summary.txt`)
   - 모든 FASTQ 파일의 분석 결과 요약

## 분석 흐름

1. FASTQ 파일 읽기
2. Indicator 시퀀스 기반 시퀀스 크롭
3. Reference 시퀀스와 정렬(alignment)
4. Edit range 내 변이 패턴 분석
5. 편집 타입 분류 및 카운팅
6. 리포트 생성

## 편집 타입 분류

- **WT**: Wild Type, 변이 없음
- **Insertion**: 염기 삽입
- **Deletion**: 염기 결실
- **Substitution**: 염기 치환
- **Complex**: 복합 편집 (여러 타입의 조합)

## 설정

`MINIMUM_FREQUENCY` 변수를 조정하여 최소 리드 카운트 임계값을 설정할 수 있습니다 (기본값: 2).

## 라이선스

Made by Chanju

## 문의

이슈나 문의사항은 GitHub Issues를 통해 제출해주세요.
