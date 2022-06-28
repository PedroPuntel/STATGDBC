Autor: Pedro Puntel (pedro.puntel@gmail.com)
Descrição: Bloco de Notas - TCC
Encoding: UTF-8

-----------
. Deadlines
-----------
	--> 04/07: 
		(STATGDBC) Implementação das correções + resultados brutos
		(Monografia) Finalizar capítulo de estatística espacial

-------
. TO-DO
-------

	-----------
	. STATGDBC:
	-----------
		--> Debugging: Motivo dos resultados com (ICS = NULL) e (grid.method = NULL)
		--> DGBClust: Utilizar matriz de distãncias original dos dados para cálculo dos índices de avaliação (silhuete e calinski)
		--> DGBClust: Questão dos testes estatísticos de hopkins/clark-evans: Mesclagem quando houver regularidade e/ou aleatoriedade espacial
		--> ESG/ASG: Verificar questão do BRKGA na comparação da probabilidade de herença (operadores de cruzamento)
		--> ESG/ASG: Verificar se a restrição de divisão das grades está funcional (principalmente para o ASG)
		--> DBSCAN: Estudar como utilizá-lo e já calcular também os seus resultados para ir comparando
	
	-------------
	. Monografia:
	-------------
	
		------------------------------------------------ Capítulo 3 ------------------------------------------------
		--> Referência Estatístca Espacial: Livro Baddley
		--> Não tem problema utilizar |M|-1 pontos no Teste de Hopkins
		--> Teste de Clark-Evans unilateral
		--> Mesclagem das células: Não queremos rejeitar a hipótese-nula de aleatoriedade espacial nos testes de Hopkins/Clark-Evans
		--> Falar brevemente sobre efeitos de borda (não considerados) e distâncias evento-evento, ponto-evento e vizinho mais próximo
		--> Descartar: Estacionariedade, Isotropia, PPP Não-Homogêneo, Monte-Carlo
		
		- Introdução:
			(OK) --> Estatística Espacial e suas subáreas/motivações: Dados de Área, Geoestatística e Processos Ponutais
			
		- Seção 3.3 - Processos Pontuais Espaciais
		    --> Mecanismo estocástico gerador dos dados, problemas de interesse, tipos de padrões de pontos
			--> Conceitos de Estacionariedade, Isotropia
			--> Definição formal + exemplos
			--> Principais processos pontuais: Processo de Poisson Homogêneo (regular e Irregular), CSR, Processo de Poisson Não-Homogêneo

		------------------------------------------------ Capítulo 5 (REVER) ------------------------------------------------
	
		- Seção 5.1 - Bases de Dados
		    --> Introdução sobre como os experimentos foram divididos, especificações computador, software e pacotes do R
			--> Convenção utilizada para com os parâmetros do BRKGA
			--> Aplicação de padronização e redução de dimensionalidade em todas as bases através do EMD
		
		- Seção 5.2.1 - Estudo de Estabilidade
			--> Tabela Resumo: Variabilidade em termos das Função Objetivo
			--> Tabela Resumo: Variabilidade em termos do número de grupos
			--> (Se possível) Tabela Resumo: Variabilidade dos tempos de execução
			--> Expressão de complexidade do algoritmo	
		
		- Seção 5.2.2 - Estudo e Comparação entre os Parâmetros
			--> Tabela Resumo: Número de aplicações do BFESGA x ESGBRKGA (Grade Simétrica, especificamente)
			--> Comparação: Qualidade das soluções x Utilizar ou Não Teste Estatístico (Etapa de Grade)
			--> Comparação: Número Médio de Grupos x Utilizar ou Não Teste Estatístico (Etapa de Grade)
			--> Comparação: Qualidade das soluções x Testes Estatísticos para Efeitos de 2° Ordem
			--> Gráfico Empate x Vitórias x Derrotas: ESG (especificando se foi BFESGA ou ESGBRKGA) x ASG
			--> Com base nos resultados, justificar escolha dos parâmetros para comparação com o DBSCAN (desempate por tempo de processamento)
			--> Avaliar a performance do EMD, comparando as projeções obtidas com o plot original das bases bidimensionais
			
		- Seção 5.2.3 - Comparação STATGDBC x DBSCAN
			--> Questão do tratamento aplicado aos objetos do tipo ruído
			--> Comparação conduzida fornecendo a ambos algoritmos as coordenadas provenientes do EMD
			--> Comparação: Número de grupos produzidos x Algoritmos
			--> Gráfico Empate x Vitórias x Derrotas: STATGDBC (especificando se foi BFESGA ou ESGBRKGA) x ASG

		------------------------------------------------ Capítulo 6 ------------------------------------------------			
		
		- Capítulo 6 - Conclusões
			-->  Considerações finais e Trabalhos futuros: Testes/Medidas estatísticas diferentes,
			Calibração de Parâmetros, Testes Estatísticos de Qualidade, Paralelização do algoritmo (disponibilidade Github),
			Comparação com outros métodos, Diferentes critérios de mesclagem, Validação mais robusta.