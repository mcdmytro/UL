{
    //common selection
    TCut cut_mva    = "mvaValue>0.2";

    if(decay==2317){
        TCut cut_m_ds17 = "m_ds17>2.305 & m_ds17<2.329";
        TCut cut_m_ds   = "m_dsi>1.955 & m_dsi<1.979";
        TCut cut_fin = cut_m_ds17&&cut_m_ds&&cut_mva;
    }

    if(decay==2460){
        TCut cut_m_ds60 = "m_ds60>2.448 & m_ds60<2.472";
        TCut cut_m_ds   = "m_dsw>1.955 & m_dsw<1.979";
        TCut cut_fin = cut_m_ds60&&cut_m_ds&&cut_mva;
    }

}
