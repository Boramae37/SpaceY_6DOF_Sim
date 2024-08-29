## Space Y 6자유도 시뮬레이션(ver 2.0) 사용법.
> 작성자 : 백현성

아래 ver 1.0과 거의 동일하지만 몇가지 변경점 있음

1. Julia 설치시에, 윈도우마켓에서 다운 받지 말고 사이트에서 설치파일 받는 것 추천.  자세한건 메뉴얼 참고

6 -> 수정; 코드 상에서만 변경하면 됨. 굳이 폴더 만들 필요 없음. 혹시 에러나면 ver1.0처럼 진행

2. 2단 추적 기능 추가

3. 각도 조정 가능하게(90도 외에도 가능) 기능 추가

4. 3d로 위치 그래프 그리는 기능 추가

5. 일정한 바람이 불때 이를 반영하는 기능 추가

> [!WARNING]
> 1. 몸체 지름은 단분리하더라도 변하지 않는것을 상정하고 설계되어있음
> 
> * 최초 실행시 15초 혹은 그 이상 걸릴 수도 있음. 이후 대략 7초 정도 걸림

----

## Space Y 6자유도 시뮬레이션(ver 1.0) 사용법.
> 작성자 : 안세환    
> 수정 : 오한준 (2024.06.12) / 백현성 (2024.06.27)

1. Julia 설치 (해당 프로그램은 1.10.4 버전으로 작성)
2. set/setup.jl을 실행하여 패키지 설치-CMD에 커맨드 입력(main.jl의 상위 디렉토리에서 julia set/setup.jl);생각보다 시간 꽤 걸림.
3. set/func.jl의 simulation() 함수 내부에서 Parameters 및 추력 정보 업데이트
4. main.jl에서 옵션 확인 (풍속, 추력(interpolation/calculate), 발사각)
5. 만약, 추력옵션 'interp' 설정 시 db/thrust.csv 를 추가해야함. 즉, 점화시간(디폴트 3초)동안 추력 데이터가 존재해야함. 양식은 초기 예시 참고
6. res에 날짜/시간 폴더를 만들고, 해당 경로로 저장 코드 변경 -> 수정
7. main.jl 실행 (vscode에서 Julia:Execute active File in REPL 로 실행할 것), 시뮬레이션 결과 저장
8. 분석 무한 반복

> [!note]
> - 세부적인 내용은 메뉴얼 참고.    
> - 풍속 데이터를 원한다면, 기상청에서 다운받아서 동서남북 잘 고려해서 고도별로 대입. (해당 부분은 이후 필요하면 추가하시길)    
> - 파일의 경로를 수정하지 않을 것.       
> - 파일 및 변수명을 임의로 바꾸지 않을 것.      
> - 최초 실행 시 약 8초의 실행시간이 걸리고 이후 1.5초 내외 소요 (그래프 저장 및 파일 저장까지)    
 




