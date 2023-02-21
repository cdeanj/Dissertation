library("dplyr")
library("ggplot2")
library("patchwork")
library("phyloseq")

AXIS_TEXT_SIZE = 9
TITLE_TEXT_SIZE = 9.5
TARGET = "s__Staphylococcus_aureus"

theme_set(theme_bw())

ps <- readRDS("../RDS/ps_shotgun_02192023.rds")
ps

psAureus <- subset_taxa(ps, Species == TARGET)

m <- psmelt(psAureus)
m <- m %>% mutate(PA = ifelse(Abundance > 0, 1, 0))

m$PA <- as.factor(m$PA)
levels(m$PA) <- c("No", "Yes")

panelA <- ggplot(m, aes(x = WeeksToInfection, fill = PA)) +
  geom_bar() +
  xlab("Weeks to infection") +
  ylab("Frequency") +
  facet_wrap(~CaseOrControl) +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = AXIS_TEXT_SIZE),
    axis.title.y = element_text(size = AXIS_TEXT_SIZE),
    plot.title = element_text(hjust = 0.5, size = TITLE_TEXT_SIZE)
  ) +
  ggtitle("Presence of Staphylococcus aureus on the teat apex")

panelB <- ggplot(m, aes(x = DaysToInfection, fill = PA)) +
  geom_bar() +
  xlab("Days to infection") +
  ylab("Frequency") +
  facet_wrap(~CaseOrControl) +
  labs(fill = "Presence") +
  theme(
    legend.position = "right",
    axis.title.x = element_text(size = AXIS_TEXT_SIZE),
    axis.title.y = element_text(size = AXIS_TEXT_SIZE)
  )

panelC <- ggplot(m, aes(x = DIM, fill = PA)) +
  geom_bar() +
  xlab("Days relative to calving") +
  ylab("Frequency") +
  facet_wrap(~CaseOrControl) +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = AXIS_TEXT_SIZE),
    axis.title.y = element_text(size = AXIS_TEXT_SIZE)
  )

figure1 <- panelA / panelB / panelC + plot_annotation(tag_levels = c("A"))
figure1

ggsave("../Figures/Figures/Figure1.png", figure1, width = 6.0, height = 9.0)
