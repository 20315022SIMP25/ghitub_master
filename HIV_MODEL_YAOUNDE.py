/**
* Name: SuiviInfecions
* Based on the internal empty template. 
* Author: Pasco
* Tags: 
*/


model SuiviInfecions

global{
	/*PARAMETRAGE DU MODELE */
	
	/* Variables du modèle */
	
	int nb_s1	; // Le nombre intiale des susceptibles moins exposés (PG)
	int nb_s2	; // Le nombre intiale des population clés ( TS)
	int nb_i	; // Le nombre initiale des infectés infectieux (I)
	int nb_c	; // Le nombre initiale des infectés Chroniques (C)
	int nb_a	; // Le nombre initiale des Infectés aiguës ou malade du sida (A)
	
	/* Parmètres du modèle */
	
	float lambda_1	; // Le recrutement dans la population générale
	float lambda_2	; // Le recrutement dans la population clé
	float mu		; // Le taux de mortalité naturelle
	float beta_1	; // Le taux d'infection du VIH dans la population générale
	float beta_2	; // Le taux d'infection du VIH dans la population clés
	float gamma		; // Le taux de d'amélioration d'état de santé des infectés I pour le compartiment des Chroniques
	float sigma		; // Le taux de migrantion des populations générale vers la population clé 
	float phi		; // Le taux de défaut traitement pour les individus C
	float alpha		; // le taux de défaut traitement pour les individus I
	float delta		; // Le taux de traitement des malades du Sida pour s'améliorer 
	float eta		; // Le taux de mortalité liée au Sida
	float epsilon	; // Le taux de retour des populations clé vers la population générale
	float rho_1		; // La sensibilisation sur le VIH/Sida
	float rho_2		; // L'usage des préservatif
	float rho_3		; // Le Test et dépistage
	
	int N <- nb_s1 + nb_s2 + nb_i + nb_c + nb_a;// Dynamique du modèle au début de l'épidémie
	
	float hKR4 <- 0.7;
	
	init{
		// CREATION DES COMPARTIMENTS DU MODELE
		
			// Création de la classe de population moins à risque dite population générale
		create S1_agt{
			S1size 			<- 	float(nb_s1)	;
			self.beta_1 	<- 	myself.beta_1	;
			self.lambda_1 	<- 	myself.lambda_1	;
			self.mu 		<- 	myself.mu		;
			self.sigma 		<- 	myself.sigma	;
			self.epsilon 	<- 	myself.epsilon	;
		}
			write 'Basic Reproduction Number (R0): ' + string((1-rho_1)*(1-rho_2)*beta_2 * (mu + phi)*(mu + eta + delta) /
			((mu + phi)*((mu + gamma*rho_3 + alpha)*(mu + eta + delta)+ alpha * delta) + (mu + eta + delta)*phi*(1-rho_1)*alpha))+'\n';
			
			// Création de la calsse de population plus à risque dite population clé
		
		create S2_agt{
			S2size 			<- 	float(nb_s2)	;
			self.beta_2 	<- 	myself.beta_2	;
			self.lambda_2 	<- 	myself.lambda_2	;
			self.mu 		<- 	myself.mu		;
			self.sigma 		<- 	myself.sigma	;
			self.epsilon 	<- 	myself.epsilon	;
			self.rho_1		<- 	myself.rho_1	;
			self.rho_2		<- 	myself.rho_2	;
		}
			// Création de la classe des infectés infectieux
		
		create I_agt {
			Isize 		<- float(nb_i)			;
			self.beta_1 <- myself.beta_1		;
			self.beta_2 <- myself.beta_2		;
			self.gamma  <- myself.gamma			;
			self.alpha	<- myself.alpha			;
			self.mu		<- myself.mu			;
			self.phi	<- myself.phi			;
			self.delta	<- myself.delta			;
			self.rho_1	<- 	myself.rho_1		;
			self.rho_2	<- 	myself.rho_2		;
			self.rho_3 	<- myself.rho_3			;
		}
			// Création de la classe des infectés chroniques
		create C_agt {
			Csize		<- float(nb_c)   ;
			self.gamma	<- myself.gamma  ;
			self.phi	<- myself.phi    ;
			self.mu		<- myself.mu     ;
			self.rho_3 	<- myself.rho_3	 ;
		}
			// Création de la classe des infectés aigues, ou symptotique encore malade du sida
		create A_agt {
			Asize		<- float(nb_a)	 	;
			self.alpha	<- myself.alpha		;
			self.delta	<- myself.delta		;
			self.eta	<- myself.eta		;
			self.mu		<- myself.mu		;

		}
			// Création du modèle SICA du système 
		create SICA_agt {
			
			self.S1 <- float(nb_s1)	;
			self.S2 <- float(nb_s2)	;
			self.I 	<- float(nb_i)	;
			self.C 	<- float(nb_c)	;
			self.A 	<- float(nb_a)	;
			
			self.beta_1 <- myself.beta_1		;
			self.beta_2 <- myself.beta_2		;
			self.gamma  <- myself.gamma			;
			self.alpha	<- myself.alpha			;
			self.mu		<- myself.mu			;
			self.phi	<- myself.phi			;
			self.delta	<- myself.delta			;
			self.lambda_1 	<- 	myself.lambda_1	;
			self.lambda_2 	<- 	myself.lambda_2	;
			self.sigma 		<- 	myself.sigma	;
			self.epsilon 	<- 	myself.epsilon	;
			self.eta		<- 	myself.eta		;
			self.rho_1 		<- 	myself.rho_1	;
			self.rho_2 		<- 	myself.rho_2	;
			self.rho_3 		<- 	myself.rho_3	;
		}
	}
}
	// La spécification des classes du modèle
	
	species S1_agt {
		
		float t		 	;
		float S1size 	; // Le nombre de susceptibles moins exposés 
		
		float beta_1	;
		float lambda_1	;
		float mu		;
		float sigma		;
		float epsilon	;
		
		// Dynamique au niveau de la population générale
		equation evol simultaneously: [  ( S2_agt ) ,  ( I_agt ), ( C_agt ), ( A_agt ) ]{
			diff ( first ( S1_agt ) . S1size , t ) = 
			(( lambda_1 - beta_1 * first ( S1_agt ) . S1size * first ( I_agt ) . Isize / N )
			 - mu * first (S1_agt) . S1size - sigma * first (S1_agt).S1size + epsilon * first (S2_agt).S2size) ;
		}
		
		reflex solving {solve evol method: "rk4" step_size: 0.001 ;} 
		
	}
	
	species S2_agt {
		
		float t			;
		float S2size 	; // Le nombre de susceptibles à risque élevé
		
		float beta_2	;
		float lambda_2	;
		float mu		;
		float sigma		;
		float epsilon	;
		float rho_1		;
		float rho_2		;
		
		
		// Dynamique au niveau de la population clé
		
		equation evol simultaneously: [  ( S1_agt ) ,  ( I_agt ), ( C_agt ), ( A_agt ) ]{
			diff ( first ( S2_agt ) . S2size , t ) = 
			(( lambda_2 - (1 - rho_1)*(1 - rho_2) * beta_2 * first ( S2_agt ) . S2size * first ( I_agt ) . Isize / N )
			 - mu * first (S2_agt) . S2size + sigma * first ( S1_agt ).S1size - epsilon * first ( S2_agt ).S2size) ;
		}
		
	}
	
	
	species I_agt {
		
		float t			;
		float Isize 	; // Le nombre des infectieux
		
		float beta_1	;
		float beta_2	;
		float gamma		;
		float alpha		;
		float mu		;
		float phi		;
		float delta		;
		float rho_1		;
		float rho_2		;
		float rho_3		;
		
		// Evolution des infectés du modèle
		
		equation evol simultaneously: [  ( S1_agt ) , ( S2_agt ), ( C_agt ), ( A_agt ) ]{
			
			diff ( first ( I_agt ). Isize , t ) =
			(((beta_1 * first ( S1_agt ) . S1size + (1 - rho_1)*(1 - rho_2) * beta_2 * first ( S2_agt ) . S2size) * first ( I_agt ) . Isize /N ) 
				- (mu + gamma * rho_3 + alpha) * first ( I_agt ) . Isize 
				+ phi * first ( C_agt ) . Csize + delta* first ( A_agt ). Asize );
		}
		
	}
	
	species C_agt {
		
		float t		 ;
		float Csize  ; // Le nombre des infectés chroniques
		
		float gamma  ;
		float phi    ;
		float mu     ;
		float rho_3	;
		
		equation evol simultaneously: [  ( S1_agt ) , ( S2_agt ), ( I_agt ), ( A_agt ) ]{
			
			diff( first ( C_agt ) .Csize, t )= 
			(rho_3 * gamma * first ( I_agt ) . Isize - ( mu + phi) * first ( C_agt ). Csize);
			
		}
	}
	
	species A_agt {
		
		float t		;
		float Asize ; 
		
		float alpha	;
		float delta	;
		float eta	;
		float mu	;
		
		
		
		equation evol simultaneously: [  ( S1_agt ) , ( S2_agt ), ( I_agt ), ( C_agt ) ]{
			
			diff( first ( A_agt ) . Asize, t )= 
			(alpha * first ( I_agt ) . Isize - ( mu + eta + delta ) * first ( A_agt ). Asize);
			
		}
		
	}
	
	// Le modèle épidémiologique mathématique SICA
	
	species SICA_agt {
		
		float t;
		
		float S1	;
		float S2	;
		float I		;
		float C		;
		float A		;
			
		float beta_1	;
		float beta_2	;
		float gamma		;
		float alpha		;
		float mu		;
		float phi		;
		float delta		;
		float lambda_1	;
		float lambda_2	;
		float sigma		;
		float epsilon	;
		float eta		;
		float rho_1		;
		float rho_2		;
		float rho_3		;
		
		// La dynamique du modèle
		
		equation SICA {
			
			diff ( S1 , t ) = (lambda_1 - beta_1 * S1 * I / N - mu * S1 - sigma * S1 + epsilon * S2) ;
		 	diff ( S2 , t ) = (lambda_2 - (1-rho_1) * (1-rho_2) * beta_2 * S2 * I / N - mu * S2 + sigma * S1 - epsilon * S2) ;
		 	diff ( I , t )  = ((beta_1 * S1 + (1-rho_1) * (1-rho_2) * beta_2 * S2) * I / N - (mu + gamma * rho_3 + alpha) * I + phi * C + delta * A );
		 	diff ( C , t )  = (gamma * rho_3 * I - ( mu + phi ) * C);
		 	diff ( A , t )  = (alpha * I - ( mu + eta + delta ) * A);
		}
		reflex solving {solve SICA method: "rk4" step_size: 0.001 ;}
	}
	
	experiment d_vih_yde type: gui{
		
		parameter 'Population générale' type: int var: nb_s1 <- 343826 category: "Initial population";
		parameter 'Population clé' type: int var: nb_s2 <- 2587 category: "Initial population";
		parameter 'Infectés' type: int var: nb_i <- 21 category: "Initial population";
		parameter 'Chroniques' type: int var: nb_s1 <- 0 category: "Initial population";
		parameter 'Aiguës' type: int var: nb_s1 <- 0 category: "Initial population";
		
		
		parameter 'Beta_1 (S1->I)'  	type: float var: beta_1 <- 0.19   		category: "Parameters";
		parameter 'Beta_2 (S2->I)'  	type: float var: beta_2 <- 2.66			category: "Parameters";
		parameter 'Gamma (I->C)'  		type: float var: gamma <- 1.0			category: "Parameters";
		parameter 'Alpha (I->A)' 		type: float var: alpha <- 0.1			category: "Parameters";
		parameter 'Mu (S,I,C,A)'  		type: float var: mu <- 0.008 			category: "Parameters";
		parameter 'Phi (C->I)'  		type: float var: phi <- 0.420 			category: "Parameters";
		parameter 'Delta (A->I)'  		type: float var: delta <- 0.034			category: "Parameters";
		parameter 'Lambda_1 (S1)'  	 	type: float var: lambda_1 <- 12378.0    category: "Parameters";
		parameter 'Lambda_2 (S2)'   	type: float var: lambda_2 <- 105.0 		category: "Parameters";
		parameter 'Sigma (S1->S2)'  	type: float var: sigma <- 0.08 			category: "Parameters";
		parameter 'Epsilon(S2->S1)' 	type: float var: epsilon <- 0.1 		category: "Parameters";
		parameter 'Eta (A)'  			type: float var: eta <- 1.0 			category: "Parameters";
		parameter 'rho_1 (S2,I)'  		type: float var: rho_1 <- 0.75 			category: "Parameters";
		parameter 'rho_2 (S2,I,C,A)'  	type: float var: rho_2 <- 0.7 			category: "Parameters";
		parameter 'rho_3 (I,C)'  		type: float var: rho_3 <- 0.35			category: "Parameters";
		
		output {
		layout #split;
		
		display "Infections" axes: false{
				
				chart 'Ctrl_Inf.Ydé ' 	 type: series background: #white {
				//data 'Susceptibles (-) ' value: first( SICA_agt ).S1 color: rgb(#green)marker:false;
				//data 'Susceptibles (+) ' value: first( SICA_agt ).S2 color: rgb(#limegreen)marker:false;
				data 'Infectieux' 		 value: first( SICA_agt ).I  color: rgb(#red)marker:false;
				//data 'Chroniques'        value: first( SICA_agt ).C  color: rgb(#maroon)marker:false;
				//data 'Aiguës'        	 value: first( SICA_agt ).A  color: rgb(#royalblue)marker:false;
				}
			}
		}
	}
