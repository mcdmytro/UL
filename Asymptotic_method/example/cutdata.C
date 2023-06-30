{
//common selections
TCut cut01="abs(lp1_dr)<0.5&&abs(lm1_dr)<0.5";
TCut cut02="abs(lp1_dz)<2.0&&abs(lm1_dz)<2.0";
TCut cut03="lp1_pt>0.1&&lm1_pt>0.1";
TCut cut04="lp1_pid<0.95&&lm1_pid<0.95";
TCut cut05="ntrk>4&&hadb_flg>0";
//TCut cut05="ntrk>2";
TCut cut06=cut01&&cut02&&cut03&&cut04&&cut05;

TCut cut30="ntrk>2";
TCut cut31=cut01&&cut02&&cut03&&cut04&&cut30;
//pip
TCut cut07="abs(pip_dr)<0.5&&abs(pip_dz)<2.0&&pip_pt>0.1";
TCut cut08="pip_eid<0.95&&pip_muid<0.95";

TCut cut09="abs(mll-3.097)<0.03";

TCut cut10="mll>2.97&&mll<3.03";
TCut cut11="mll>3.17&&mll<3.23";
TCut cut12=cut10||cut11;

TCut cut13="lp1_eid>0.95&&lm1_eid>0.95";
TCut cut14="lp1_eid>0.05&&lm1_eid>0.05";
TCut cut15="lp1_muid>0.95||lm1_muid>0.95";
TCut cut16="lp1_muid>0.05&&lm1_muid>0.05";
TCut cut17=cut13&&cut14;
TCut cut18=cut15&&cut16;
TCut cut19=cut17||cut18;

TCut cut20="mgfit01>3.8&&mgfit01<4";
TCut cut21="refit01>3.857&&refit01<3.917";

//sidebands

TCut cut22="mgfit01>3.&&mgfit01<5";
}
