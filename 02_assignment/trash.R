We now replicate table 2
```{r}
list_2 <- dat2 %>%
  select(COUNTRY, geometry, NAME_1, NAME_2) %>%
  mutate(country = ifelse(COUNTRY == "Brazil", "BRA", COUNTRY)) %>%
  rename("muni" = "NAME_2") %>%
  mutate(state = ifelse(NAME_1 == "Rio Grande do Sul", "RS", NAME_1)) 

literacy_Arg_Bra_Par_2 <- literacy_Arg_Bra_Par %>%
  left_join(list_2, by = c("muni",  "state", "country"))

mod1 <- lm(illiteracy ~ (lati) + (longi) + distmiss + state, data = literacy_Arg_Bra_Par_2)

mod2 <- lm(illiteracy ~ (lati) + (longi) + distmiss + state + coast + river + slope + rugg + alti + tempe + preci + area, data = literacy_Arg_Bra_Par_2)


bra <- literacy_Arg_Bra_Par_2 %>% 
  filter(country == "BRA")

mod3 <- lm(illiteracy ~ (lati) + (longi) + distmiss + mesorregi, data = bra)

mod4 <- lm(illiteracy ~ (lati) + (longi) + distmiss + mesorregi + coast + river + slope + rugg + alti + tempe + preci + area, data = bra)


arg <- literacy_Arg_Bra_Par_2 %>% 
  filter(country == "Argentina")

mod5 <- lm(illiteracy ~ (lati) + (longi) + distmiss, data = arg)

mod6 <- lm(illiteracy ~ (lati) + (longi) + distmiss + coast + river + slope + rugg + alti + tempe + preci + area, data = arg)


par <- literacy_Arg_Bra_Par_2 %>% 
  filter(country == "Paraguay")

mod7 <- lm(illiteracy ~ (lati) + (longi) + distmiss + state, data = par)

mod8 <- lm(illiteracy ~ (lati) + (longi) + distmiss + state + coast + river + slope + rugg + alti + tempe + preci + area, data = par)

#stargazer(mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8)

```
